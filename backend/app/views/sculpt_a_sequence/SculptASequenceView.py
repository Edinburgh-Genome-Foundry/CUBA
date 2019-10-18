"""Bla."""

from rest_framework import serializers
from ..serializers import FileSerializer
from dnachisel import (DnaOptimizationProblem, sequence_to_biopython_record,
                       annotate_record)


from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data)


class LabeledFeatureSerializer(serializers.Serializer):
    strand = serializers.IntegerField()
    start = serializers.IntegerField()
    end = serializers.IntegerField()
    label = serializers.CharField()


class serializer_class(serializers.Serializer):
    file = FileSerializer()
    editFeatures = serializers.BooleanField()
    editedFeatures = serializers.DictField(child=LabeledFeatureSerializer())
    sequence = serializers.CharField()

class worker_class(AsyncWorker):

    def work(self):

        data = self.data
        self.logger(message='Initializing...')

        if data.editFeatures:
            record = sequence_to_biopython_record(data.sequence.upper())
            for feature in sorted(data.editedFeatures.values(),
                                  key=lambda f: (f.start, f.end)):
                annotate_record(record, feature_type="misc_feature",
                                location=(feature.start, feature.end),
                                label=feature.label)
        else:
            record = records_from_data_files([data.file])[0]
        problem = DnaOptimizationProblem.from_record(record, logger=self.logger)
        problem.optimization_stagnation_tolerance = 30
        success, summary, zip_data = problem.optimize_with_report(
            target="@memory", project_name=record.id)
        return {
          'zip_file': {
              'data': data_to_html_data(zip_data, 'zip'),
              'name': 'optimization_report.zip',
              'mimetype': 'application/zip'
          },
          'success': success,
          'summary': summary
        }

class SculptASequenceView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
