"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import records_from_data_files, data_to_html_data
from ..serializers import FileSerializer
from genedom import batch_domestication, BUILTIN_DOMESTICATORS
import flametree


class serializer_class(serializers.Serializer):
    """Serializer."""
    parts = serializers.ListField(child=FileSerializer())
    standard = serializers.CharField()
    standard_name = serializers.CharField(allow_null=True, allow_blank=True)
    standard_definition = FileSerializer(allow_null=True, required=False)

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        records = records_from_data_files(data.parts)

        if data.standard == 'EMMA':
            def domesticator(record):
                """Find which domesticator to use for the given part.
                Here we use the fact that all records will have an
                ID of the form ``position_partname``. We extract the
                position and return the corresponding domesticator.
                """
                position = record.id.split("_")[0]
                return BUILTIN_DOMESTICATORS.EMMA[position]
        else:
            raise ValueError("Support for %s not implemented" % data.standard)
        nfails, zip_data = batch_domestication(
            records, domesticator, '@memory', allow_edits=True)
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
