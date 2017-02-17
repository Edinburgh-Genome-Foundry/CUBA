"""Bla."""
from rest_framework import serializers
from .base import AsyncWorker, StartJobView
import requests
import time
from website.settings import REPORTING_HOST
from django.urls import reverse

class serializer_class(serializers.Serializer):
    name = serializers.CharField(label="Name")
    legs = serializers.IntegerField(label="Number of legs")
    eyes = serializers.IntegerField(label="Number of eyes")

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        print("i am here AH AH AH")
        print(self.data)

        self.set_progress_message("Job (fakely) in progress")
        text = (
            "Congratulations! your new pet has been generated. It is called "
            "{name}, has {legs} legs and {eyes} eyes."
        ).format(**self.data)
        time.sleep(1)

        self.set_progress_message("(Fake) report generation")
        request_data = {
            "text": text,
            "generate_report": self.generate_report,
            "generate_preview": self.generate_preview,
        }
        reports_uri = "http://%s%s" % (REPORTING_HOST + ':8000',
                                       reverse("reports"))
        response = requests.post(reports_uri, request_data)
        print(reports_uri)
        print(response)
        time.sleep(1)

        return response.json()


class ScenarioOneView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
