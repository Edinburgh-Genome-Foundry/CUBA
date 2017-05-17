"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from collections import OrderedDict
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import string_to_record, matplotlib_figure_to_svg_base64_data
from bandwitch import (IdealDigestionsProblem, SeparatingDigestionsProblem,
                       NoSolutionError, LADDERS)


digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()
    circularity = serializers.BooleanField()

class ComaSeparatedList(serializers.CharField):

    def to_representation(self, value):
        return [e.strip() for e in value.split(",")]

class serializer_class(serializers.Serializer):
    goal = serializers.CharField()
    bandsRange = serializers.ListField(child=serializers.IntegerField())
    files = serializers.ListField(child=FileSerializer())
    ladder = serializers.CharField()
    bandsPrecision = serializers.IntegerField()
    maxEnzymes = serializers.IntegerField()
    possibleEnzymes = ComaSeparatedList()

class worker_class(AsyncWorker):

    def work(self):

        self.set_progress_message('Exploring possible digestions...')

        data = self.data
        ladder = LADDERS[data.ladder]
        band_min, band_max = data.bandsRange

        sequences = OrderedDict()
        for f in data.files:
            content = f.content.split("base64,")[1]
            content = b64decode(content).decode("utf-8")
            record, fmt = string_to_record(content)
            sequences[f.name] = str(record.seq)

        if (data.goal == 'ideal'):
            self.set_progress_message("Finding good digestions...")
            class MyIdealDigestionsProblem(IdealDigestionsProblem):
                def migration_pattern_is_ideal(self, migration):
                    min_m = 0.9 * self.migration_min + 0.1 * self.migration_max
                    max_m = 0.1 * self.migration_min + 0.9 * self.migration_max
                    bands_in_central_zone = [band for band in migration
                                             if min_m <= band <= max_m]
                    return band_min <= len(bands_in_central_zone) <= band_max

            # DEFINE THE SEQUENCES AND THE ENZYME SET

            problem = MyIdealDigestionsProblem(
                 sequences, data.possibleEnzymes, linear=False, ladder=ladder,
                 max_enzymes_per_digestion=data.maxEnzymes
            )
        else:
            print (0.01 * data.bandsPrecision)
            problem = SeparatingDigestionsProblem(
                 sequences, data.possibleEnzymes, linear=False, ladder=ladder,
                 relative_error=0.01 * data.bandsPrecision,
                 max_enzymes_per_digestion=data.maxEnzymes
            )
        try:
            self.set_progress_message('Selecting digestions...')
            selected_digestions = problem.select_digestions()
            axes = problem.plot_digestions(
                selected_digestions,
                patterns_props={'label_fontdict': {'rotation': 35}}
            )
            figure_data = matplotlib_figure_to_svg_base64_data(
                axes[0].figure, bbox_inches='tight'
            )
            return {
              'figure_data': figure_data,
              'digestions': selected_digestions
            }
        except NoSolutionError:
            return {
              'digestions': []
            }

class SelectDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
