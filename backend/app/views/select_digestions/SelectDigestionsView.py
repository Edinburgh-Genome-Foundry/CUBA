"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from collections import OrderedDict
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files,
                     matplotlib_figure_to_svg_base64_data)
from bandwitch import (IdealDigestionsProblem, SeparatingDigestionsProblem,
                       LADDERS)


digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class ComaSeparatedList(serializers.CharField):
    def to_representation(self, value):
        return [e.strip() for e in value.split(",")]

class serializer_class(serializers.Serializer):
    goal = serializers.CharField()
    bandsRange = serializers.ListField(child=serializers.IntegerField())
    files = serializers.ListField(child=FileSerializer())
    ladder = serializers.CharField()
    maxDigestions = serializers.IntegerField()
    maxEnzymes = serializers.IntegerField()
    circularSequences = serializers.BooleanField()
    possibleEnzymes = ComaSeparatedList()

class worker_class(AsyncWorker):

    def work(self):

        self.set_progress_message('Exploring possible digestions...')

        data = self.data
        ladder = LADDERS[data.ladder]
        enzymes = data.possibleEnzymes

        sequences = OrderedDict([
            (record.id, str(record.seq))
            for record in records_from_data_files(data.files)
        ])
        # print (sequences)
        self.set_progress_message("Initializing...")
        if (data.goal == 'ideal'):

            mini, maxi = data.bandsRange

            problem = IdealDigestionsProblem(
                sequences=sequences, enzymes=enzymes, ladder=ladder,
                linear=not data.circularSequences,
                min_bands=mini, max_bands=maxi,
                max_enzymes_per_digestion=data.maxEnzymes
            )
        else:
            print(enzymes)
            problem = SeparatingDigestionsProblem(
                 sequences=sequences, enzymes=enzymes, ladder=ladder,
                 linear=not data.circularSequences,
                 max_enzymes_per_digestion=data.maxEnzymes
            )
        self.set_progress_message('Selecting digestions...')
        score, selected_digestions = problem.select_digestions(
            max_digestions=data.maxDigestions, search='full')
        axes = problem.plot_digestions(
            selected_digestions,
            patterns_props={'label_fontdict': {'rotation': 35}}
        )
        figure_data = matplotlib_figure_to_svg_base64_data(
            axes[0].figure, bbox_inches='tight'
        )
        return {
          'figure_data': figure_data,
          'digestions': selected_digestions,
          'score': score,
        }

class SelectDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
