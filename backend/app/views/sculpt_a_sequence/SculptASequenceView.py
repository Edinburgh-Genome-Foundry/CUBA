"""Bla."""

from rest_framework import serializers
from dnachisel.reports import optimization_with_report
from dnachisel import DnaOptimizationProblem

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_file,
                     zip_data_to_html_data)

class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class serializer_class(serializers.Serializer):
    file = FileSerializer()

class worker_class(AsyncWorker):

    def work(self):

        data = self.data
        records, fmt = records_from_data_file(data.file)
        record = records[0]
        problem = DnaOptimizationProblem.from_record(record)
        problem.progress_logger = self.update_progress_data
        success, summary, zip_data = optimization_with_report(
            target="@memory", problem=problem, project_name=record.id)
        return {
          'zip_file': {
              'data': zip_data_to_html_data(zip_data),
              'name': 'optimization_report.zip',
              'mimetype': 'application/zip'
          },
          'success': success,
          'summary': summary
        }

class SculptASequenceView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
