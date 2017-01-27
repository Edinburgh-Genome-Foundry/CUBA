"""Bla."""
from rest_framework import serializers
from .base import AsyncWorker, StartJobView
import requests
import time
from django.urls import reverse

class serializer_class(serializers.Serializer):
    name = serializers.CharField(label="Name")
    legs = serializers.IntegerField(label="Number of legs")
    eyes = serializers.IntegerField(label="Number of eyes")

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        self.set_progress_message("Job (fakely) pending")
        time.sleep(3)

        self.set_progress_message("Job (fakely) in progress")
        text = ("Congratulations, your new pet has been generated. It's called"
                "%{name}, has %{legs} legs and %{eyes} eyes.") % self.data
        time.sleep(3)

        self.set_progress_message("(Fake) report generation")
        request_data = {
            "text": text,
            "generate_report": self.generate_report,
            "generate_preview": self.generate_preview,
        }
        reports_uri = "http://%s%s" % (self.domain_name, reverse("reports"))
        response = requests.post(reports_uri, request_data)
        time.sleep(3)

        self.set_progress_message("Done ! Results are (fakely) on their way !")
        time.sleep(2)

        return response




class ScenarioTwoView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
