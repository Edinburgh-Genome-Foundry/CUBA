import json

from rest_framework_tracking.mixins import LoggingMixin

from rest_framework import status, serializers
from rest_framework.response import Response
from rest_framework.views import APIView
from rest_framework.renderers import JSONRenderer

from django.urls import reverse

import django_rq
import rq

import dnaweaver as dw
import requests

class ObjectDict(dict):

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        dict.__getattribute__(self, key)

    def __setattr__(self, key, value):
        self[key]=value
        self.__dict__[key] = value

    @staticmethod
    def from_dict(d):
        obj = ObjectDict({
            key: ([
                ObjectDict.from_dict(e)
                if isinstance(e, dict) else e
                for e in value
            ]
                if isinstance(value, (list, tuple))
                else (ObjectDict.from_dict(value)
                      if isinstance(value, dict)
                      else value))
            for (key, value) in d.items()
        })
        for key, value in obj.items():
            sanitized_key = key.replace(" ", "_").replace(".", "_")
            obj.__dict__[sanitized_key] = value
        return obj


class SerializerView(APIView):
    def serialize(self, request):
        if hasattr(self, 'serializer_class'):
            serializer = self.serializer_class(data=request.data)
            if serializer.is_valid():
                return ObjectDict.from_dict(serializer.data)
            else:
                return Response(serializer.errors,
                                status=status.HTTP_400_BAD_REQUEST)
        else:
            return ObjectDict.from_dict(request.data)




class PollJobView(SerializerView):
    renderer_classes = (JSONRenderer,)
    class serializer_class(serializers.Serializer):
        job_id = serializers.CharField(label="Job ID")

    def post(self, request, format=None):
        """A view to report the progress to the user."""
        data = self.serialize(request)
        # job_id = request.POST['job_id']
        job = django_rq.get_queue("default").fetch_job(data.job_id)
        if job is None:
            return Response(dict(success=False, error="Unknown job ID."))
        return Response(
            dict(
                success=True,
                status=job.get_status(),
                error='',
                progress=dict(
                  data=job.meta.get('progress_data', None),
                  message=job.meta.get('progress_message', None)
                ),
                result=job.result
            ),
            status=status.HTTP_200_OK
        )

class StartJobView(SerializerView, LoggingMixin):
    renderer_classes = (JSONRenderer, )
    def post(self, request, format=None):
        """ Some description for posts"""
        data = self.serialize(request)
        if not isinstance(data, dict):
            return data
        data["domain_name"] = request.META['HTTP_HOST']
        job = django_rq.get_queue("default").enqueue(
            self.worker_class.run, data)
        return Response({"job_id": job.id}, status=status.HTTP_200_OK)

class JobResult:

    def __init__(self, json=None, preview_html=None, preview_js=None,
                 file_data=None, file_name=None, file_mimetype=None,
                 error=None):
        self.error = error
        self.json = json
        self.preview_html = preview_html
        self.preview_js = preview_js
        self.file_data = file_data
        self.file_mimetype = file_mimetype
        self.file_name = file_name

    def as_json(self):
        return dict(
            error=self.error,
            json=self.json,
            preview=dict(html=self.preview_html, js=self.preview_js),
            file=dict(data=self.file_data, name=self.file_name,
                      mimetype=self.file_mimetype)
        )

class AsyncWorker:

    def __init__(self, data, job):

        self.job = job
        self.data = data
        self.domain_name = self.data.pop("domain_name", None)

    def set_progress_message(self, message):
        self.job.meta['progress_message'] = message
        self.job.save()

    def set_progress_data(self, data):
        self.job.meta['progress_data'] = data
        self.job.save()

    @classmethod
    def run(cls, data):
        # worker = django_rq.get_worker('default')
        job = rq.get_current_job()
        agent = cls(data, job)
        result = agent.work()
        if isinstance(result, JobResult):
            result = result.as_json()
        return result
