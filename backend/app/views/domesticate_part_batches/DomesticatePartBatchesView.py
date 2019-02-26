"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data,
                     spreadsheet_file_to_dataframe)
from ..serializers import FileSerializer
from genedom import (batch_domestication, BUILTIN_STANDARDS,
                     GoldenGateDomesticator)
import flametree


class serializer_class(serializers.Serializer):
    """Serializer."""
    parts = serializers.ListField(child=FileSerializer())
    standard = serializers.CharField()
    standard_name = serializers.CharField(allow_null=True, allow_blank=True)
    standard_definition = FileSerializer(allow_null=True, required=False)
    allow_edits = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        records = records_from_data_files(data.parts)
        self.logger(message="Now domesticating the parts")

        if data.standard == 'EMMA':
            standard = BUILTIN_STANDARDS.EMMA
        elif data.standard == 'custom':
            dataframe = spreadsheet_file_to_dataframe(data.standard_definition)
            standard = GoldenGateDomesticator.standard_from_spreadsheet(
                dataframe=dataframe, name_prefix=data.standard_name
            )
        nfails, zip_data = batch_domestication(
            records, '@memory', standard=standard, logger=self.logger,
            allow_edits=data.allow_edits)
        return dict(
          file=dict(data=data_to_html_data(zip_data, 'zip'),
                    name='domestication_report.zip',
                    mimetype='application/zip'),
          nfails=nfails,
          success=nfails == 0
        )




class DomesticatePartBatchesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
