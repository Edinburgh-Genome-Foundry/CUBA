"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import data_to_html_data, records_from_data_file
from bandwitch import LADDERS
import bandwagon
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
    makeReport = serializers.BooleanField()
    showBandsSizes = serializers.BooleanField()

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data
        records = [
            records_from_data_file(f)[0][0]
            for f in data.files
        ]
        for f, record in zip(data.files, records):
            record.linear = not f.circularity
            record.id = f.name
        ladder = bandwagon.custom_ladder(None, LADDERS[data.ladder].bands)
        digestions = data.digestions

        self.logger(message="Generating report...")
        preview, report = generate_report(records=records,
                                          digestions=digestions,
                                          ladder=ladder,
                                          group_by="digestions",
                                          full_report=data.makeReport,
                                          show_band_sizes=data.showBandsSizes)
        if data.makeReport:
            report = data_to_html_data(report, 'zip')
        return JobResult(
            preview_html='<img src="%s"/>' % preview,
            file_data=report,
            file_name=None if (report is None) else "digestion_report.zip",
            file_mimetype="application/zip"
        )


class PredictDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
