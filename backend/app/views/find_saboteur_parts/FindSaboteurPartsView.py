"""Bla."""

from rest_framework import serializers
from saboteurs import csv_to_groups_data, find_saboteurs, analysis_report
from ..base import AsyncWorker, StartJobView
from ..tools import data_to_html_data, file_to_filelike_object
from ..serializers import FileSerializer


class serializer_class(serializers.Serializer):
    """Serializer."""
    assemblies_data_file = FileSerializer()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        self.logger(message='analyzing results...')
        asm_data_filelike = file_to_filelike_object(data.assemblies_data_file)
        groups_data = csv_to_groups_data(asm_data_filelike)
        analysis_results = find_saboteurs(groups_data)

        self.logger(message='Writing a report...')
        report_data = analysis_report(
            analysis_results,
            "@memory",
            replacements=[
                ('groups', 'assemblies'),
                ('group', 'assembly'),
                ('member', 'part')
            ]
        )
        return {
          'report': {
              'data': data_to_html_data(report_data, 'pdf'),
              'name': 'saboteur_report.pdf',
              'mimetype': 'application/pdf'
          }
        }

class FindSaboteurPartsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
