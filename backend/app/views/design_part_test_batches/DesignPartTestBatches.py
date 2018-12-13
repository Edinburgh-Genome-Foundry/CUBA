"""Bla."""

from rest_framework import serializers
from saboteurs import (generate_combinatorial_groups, design_test_batch,
                       generate_batch_report)
from ..base import AsyncWorker, StartJobView
from ..tools import data_to_html_data, file_to_filelike_object
from ..serializers import FileSerializer


class serializer_class(serializers.Serializer):
    """Serializer."""
    input_file = FileSerializer()
    input_format = serializers.CharField()
    max_saboteurs = serializers.IntegerField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        self.logger(message='analyzing the data...')
        asm_data = file_to_filelike_object(data.input_file).read().decode()
        lines = [[e for e in l.split(',') if len(e)]
                 for l in asm_data.split('\n')
                 if len(l)]
        if data.input_format == 'combinatorial':
            combinatorial_design = {l[0]: l[1:] for l in lines[1:]}
            groups = generate_combinatorial_groups(combinatorial_design)
        else:
            groups = {l[0]: l[1:] for l in lines[1:]}
        selected_groups, error = design_test_batch(
            groups, max_saboteurs=data.max_saboteurs)
        if error:
            raise ValueError(error)
        report_data = generate_batch_report(selected_groups)
        return {
            'n_selected': len(selected_groups),
            'report': {
                'data': data_to_html_data(report_data, 'zip'),
                'name': 'test_batch.zip',
                'mimetype': 'application/zip'
            }
        }
        

class DesignPartTestBatchesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
