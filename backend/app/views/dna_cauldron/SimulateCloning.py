"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import string_to_record
from dnacauldron import full_assembly_report, autoselect_enzyme


digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()
    circularity = serializers.BooleanField()

class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    parts = serializers.ListField(child=FileSerializer())

class worker_class(AsyncWorker):

    def work(self):
        self.set_progress_message("Reading Data...")
        data = self.data

        records = []
        for f in data.parts:
            content = f.content.split("base64,")[1]
            content = b64decode(content).decode("utf-8")
            record, fmt = string_to_record(content)
            record.name = f.name
            record.linear = not f.circularity
            records.append(record)

        if data.enzyme == "Autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)

        self.set_progress_message("Generating a report, be patient.")

        nconstructs, zip_data = full_assembly_report(
            records, target='@memory', enzyme=self.data.enzyme,
            max_assemblies=40, fragments_filters='auto',
            assemblies_prefix='assembly'
        )
        zip_data = ('data:application/zip;base64,' +
                    b64encode(zip_data).decode("utf-8"))
        if nconstructs == 0:
            preview_html = 'No possible construct found, see report for more.'
        elif nconstructs == 1:
            preview_html = '1 construct was generated.'
        else:
            preview_html = "%d constructs were generated." % nconstructs

        return JobResult(
            preview_html=preview_html,
            file_data=zip_data,
            file_mimetype='application/zip',
            file_name='asm_report.zip'
        )

class SimulateCloningView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
