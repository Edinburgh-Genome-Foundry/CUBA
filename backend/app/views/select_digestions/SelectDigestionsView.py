"""Bla."""

from rest_framework import serializers
from collections import OrderedDict
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files,
                     data_to_html_data,
                     matplotlib_figure_to_svg_base64_data)
from bandwitch import (IdealDigestionsProblem, SeparatingDigestionsProblem,
                       LADDERS)
from bandwitch.plots import plot_all_constructs_cuts_maps


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
    showBandsSizes = serializers.BooleanField()
    plotCuts = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):

        self.logger(message='Exploring possible digestions...')

        data = self.data
        ladder = LADDERS[data.ladder]
        enzymes = data.possibleEnzymes
        records = records_from_data_files(data.files)
        sequences = OrderedDict([
            (record.id, str(record.seq))
            for record in records
        ])

        self.logger(message="Initializing...")

        if (data.goal == 'ideal'):
            mini, maxi = data.bandsRange
            problem = IdealDigestionsProblem(
                sequences=sequences, enzymes=enzymes, ladder=ladder,
                linear=not data.circularSequences,
                min_bands=mini, max_bands=maxi,
                max_enzymes_per_digestion=data.maxEnzymes
            )
        else:
            problem = SeparatingDigestionsProblem(
                 sequences=sequences, enzymes=enzymes, ladder=ladder,
                 linear=not data.circularSequences,
                 max_enzymes_per_digestion=data.maxEnzymes
            )
        self.logger(message='Selecting digestions...')
        score, selected_digestions = problem.select_digestions(
            max_digestions=data.maxDigestions, search='full')
        bands_props = None if not data.showBandsSizes else dict(
            label='=size',
            label_fontdict=dict(size=6)
        )
        axes = problem.plot_digestions(
            selected_digestions,
            patterns_props={'label_fontdict': {'rotation': 35}},
            bands_props=bands_props
        )
        figure_data = matplotlib_figure_to_svg_base64_data(
            axes[0].figure, bbox_inches='tight'
        )

        if data.plotCuts:
            self.logger(message="Plotting cuts maps...")
            pdf_data = plot_all_constructs_cuts_maps([
                (rec, digestion)
                for rec in records
                for digestion in selected_digestions
            ])
            pdf_file = dict(data=data_to_html_data(pdf_data, 'pdf'),
                            name='restrictions_cuts.pdf',
                            mimetype='application/pdf')
        else:
            pdf_file = None

        return {
          'figure_data': figure_data,
          'digestions': selected_digestions,
          'score': score,
          'pdf_file': pdf_file
        }

class SelectDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
