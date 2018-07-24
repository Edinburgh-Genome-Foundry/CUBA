"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files, file_to_filelike_object,
                     data_to_html_data)
from ..serializers import FileSerializer
from dnacauldron import (full_assembly_report, autoselect_enzyme,
                         full_assembly_plan_report)
from plateo import AssemblyPlan
import pandas

digestion = serializers.ListField(child=serializers.CharField())
# class SequenceFileSerializer(FileSerializer):
#     circularity = serializers.BooleanField()

class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    parts = serializers.ListField(child=FileSerializer())
    connectors = serializers.ListField(child=FileSerializer())
    include_fragments = serializers.BooleanField()
    use_assembly_plan = serializers.BooleanField()
    single_assemblies = serializers.BooleanField()
    assembly_plan = FileSerializer(allow_null=True)

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data

        records = records_from_data_files(data.parts)
        connector_records = records_from_data_files(data.connectors)
        for r in (records + connector_records):
            if not hasattr(r, 'linear'):
                r.linear = False
        if data.enzyme == "Autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)

        self.logger(message="Generating a report, be patient.")

        if data.use_assembly_plan:
            if data.assembly_plan.name.lower().endswith('.csv'):
                dataframe = pandas.DataFrame([
                    line.split(',')
                    for line in data.assembly_plan.content.split('\n')
                ])
            else:
                filelike = file_to_filelike_object(data.assembly_plan)
                dataframe = pandas.read_excel(filelike, header=None)
            assembly_plan = AssemblyPlan.from_spreadsheet(dataframe=dataframe)
            assembly_plan.parts_data = {r.id: {'record': r} for r in records}
            parts_without_data = assembly_plan.parts_without_data()
            if len(parts_without_data):
                return {
                    'success': False,
                    'unknown_parts': parts_without_data
                }
            errors, zip_data = full_assembly_plan_report(
                assembly_plan.assemblies_with_records(), target="@memory",
                enzyme=data.enzyme,
                assert_single_assemblies=data.single_assemblies,
                logger=self.logger, connector_records=connector_records,
                fail_silently=True,
                include_fragments_plots=data.include_fragments,
                include_parts_plots=data.include_fragments
            )
            infos = dict(errors=errors)

        else:

            nconstructs, zip_data = full_assembly_report(
                records,
                connector_records=connector_records,
                target='@memory',
                enzyme=data.enzyme,
                max_assemblies=40, fragments_filters='auto',
                assemblies_prefix='assembly',
                include_fragments_plots=data.include_fragments,
                include_parts_plots=data.include_fragments

            )
            # zip_data = ('data:application/zip;base64,' +
            #             b64encode(zip_data).decode("utf-8"))
            infos = dict(nconstructs=nconstructs)

        return {
             'file': {
                 'data': data_to_html_data(zip_data, 'zip'),
                 'name': 'assemblies.zip',
                 'mimetype': 'application/zip'
             },
             'success': True,
             'infos': infos
        }

class SimulateGGAssembliesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
