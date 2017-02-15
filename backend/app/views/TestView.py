"""Bla."""
from rest_framework import serializers
from rest_framework.response import Response
from .base import AsyncWorker, StartJobView
from rest_framework import status
import time

class serializer_class(serializers.Serializer):
    name = serializers.CharField(label="Name")
    legs = serializers.IntegerField(label="Number of legs")
    eyes = serializers.IntegerField(label="Number of eyes")

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        time.sleep(1)

        return Response({"content": "content"}, status.HTTP_200_OK)


class TestView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
