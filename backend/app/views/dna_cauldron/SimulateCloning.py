"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from dnacauldron import full_assembly_report


digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()
    circularity = serializers.BooleanField()

class serializer_class(serializers.Serializer):
    ladder = serializers.CharField()
    digestions = serializers.ListField(
        child=serializers.ListField(
            child=serializers.CharField()))
    files = serializers.ListField(child=FileSerializer())
    make_report = serializers.BooleanField()

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        self.set_progress_message("Reading Data...")
        print (self.data)
        return None

class SimulateCloningView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
