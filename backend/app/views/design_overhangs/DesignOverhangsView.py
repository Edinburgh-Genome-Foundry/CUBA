"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_file, zip_data_to_html_data)
from goldenhinges import OverhangsSelector
from goldenhinges.reports import write_report_for_cutting_solution


class FileSerializer(serializers.Serializer):
    """Serializer for files."""

    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    """Serializer."""

    sequence = FileSerializer(allow_null=True)
    goal = serializers.CharField()
    overhangs_differences = serializers.IntegerField()
    gc_content = serializers.ListField(child=serializers.IntegerField())
    mandatory_overhangs = serializers.ListField(child=serializers.CharField())
    forbidden_overhangs = serializers.ListField(child=serializers.CharField())
    cutting_mode = serializers.CharField()
    extremities = serializers.BooleanField()
    n_overhangs = serializers.IntegerField()
    n_fragments = serializers.IntegerField()
    auto_overhangs = serializers.BooleanField()
    left_flank_sequence = serializers.CharField(allow_blank=True)
    right_flank_sequence = serializers.CharField(allow_blank=True)
    allow_edits = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data

        selector = OverhangsSelector(
            gc_min=1.0 * data.gc_content[0] / 100,
            gc_max=1.0 * data.gc_content[1] / 100,
            differences=data.overhangs_differences,
            forbidden_overhangs=data.forbidden_overhangs,
            time_limit=2, external_overhangs=(),
            progress_logger=self.update_progress_data
        )

        if data.goal == 'sequence_decomposition':
            records, fmt = records_from_data_file(data.sequence)
            record = records[0]
            sequence = record
            self.set_progress_message("Decomposing the sequence...")
            if data.cutting_mode == 'equal':
                solution = selector.cut_sequence(
                    equal_segments=data.n_fragments, sequence=sequence,
                    max_radius=20,
                    include_extremities=data.extremities,
                    allow_edits=data.allow_edits)
            else:
                solution = selector.cut_sequence(
                    sequence=sequence,
                    include_extremities=data.extremities,
                    allow_edits=data.allow_edits)
            print("Solution", solution)
            if solution is None:
                return {
                  'success': False
                }

            zip_data = write_report_for_cutting_solution(
                solution=solution, target="@memory", sequence=sequence,
                left_flank=data.left_flank_sequence,
                right_flank=data.right_flank_sequence
            )
            return {
              'zip_file': {
                  'data': zip_data_to_html_data(zip_data),
                  'name': 'sequence_decomposition_report.zip',
                  'mimetype': 'application/zip'
              },
              'success': True
            }
        elif data.goal == 'overhangs_set':
            self.set_progress_message("Designing a set of overhangs...")

            try:
                overhangs = selector.generate_overhangs_set(
                    n_overhangs=None,
                    mandatory_overhangs=data.mandatory_overhangs,
                    n_cliques=30000
                )
            except ValueError as err:
                return {
                  'success': False,
                  'message': " ".join(err.args)
                }

            overhangs = selector.generate_overhangs_set(
                n_overhangs=None if data.auto_overhangs else data.n_overhangs,
                start_at=len(overhangs),
                mandatory_overhangs=data.mandatory_overhangs,
            )

            if overhangs is None:
                return {
                  'success': False,
                  'overhangs': []
                }
            return {
              'overhangs': overhangs,
              'success': True,
            }
        else:
            raise ValueError("Unknown goal")




class DesignOverhangsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class