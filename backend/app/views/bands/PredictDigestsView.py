"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import string_to_record, LADDERS
from .report_generator import generate_report


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
        data = self.data
        records = []
        for f in data.files:
            content = f.content.split("base64,")[1]
            content = b64decode(content).decode("utf-8")
            record, fmt = string_to_record(content)
            if record.id == "<unknown id>":
                record.id = f.name
                record.linear = not f.circularity
            records.append(record)
        ladder = LADDERS[data.ladder]
        digestions = data.digestions

        self.set_progress_message("Generating report...")
        preview, report = generate_report(records=records,
                                          digestions=digestions,
                                          ladder=ladder,
                                          group_by="digestions",
                                          full_report=data.make_report)
        if data.make_report:
            report = ('data:application/zip;base64,' +
                      b64encode(report).decode("utf-8"))
        return JobResult(
            preview_html='<img src="%s"/>' % preview,
            file_data=report,
            file_name=None if (report is None) else "digestion_report.zip",
            file_mimetype="application/zip"
        )


class PredictDigestsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
