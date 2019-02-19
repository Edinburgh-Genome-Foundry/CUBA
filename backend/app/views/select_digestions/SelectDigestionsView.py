"""Bla."""

from rest_framework import serializers
from collections import OrderedDict
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files,
                     data_to_html_data,
                     matplotlib_figure_to_svg_base64_data)
from bandwitch import (IdealDigestionsProblem, SeparatingDigestionsProblem,
                       LADDERS)
import bandwagon
import flametree
from bandwagon import plot_records_digestions
from io import BytesIO


digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class ComaSeparatedList(serializers.CharField):
    def to_representation(self, value):
        return [e.strip() for e in value.split(",")]

class serializer_class(serializers.Serializer):
    goal = serializers.CharField()
    bands_range = serializers.ListField(child=serializers.IntegerField())
    files = serializers.ListField(child=FileSerializer())
    ladder = serializers.CharField()
    max_digestions = serializers.IntegerField()
    max_enzymes = serializers.IntegerField()
    circular_sequences = serializers.BooleanField()
    possible_enzymes = ComaSeparatedList()
    show_bands_sizes = serializers.BooleanField()
    plot_cuts = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):

        self.logger(message='Exploring possible digestions...')

        data = self.data
        ladder = LADDERS[data.ladder]
        enzymes = data.possible_enzymes
        records = records_from_data_files(data.files)
        sequences = OrderedDict([
            (record.id, str(record.seq))
            for record in records
        ])

        self.logger(message="Initializing...")

        if (data.goal == 'ideal'):
            mini, maxi = data.bands_range
            problem = IdealDigestionsProblem(
                sequences=sequences, enzymes=enzymes, ladder=ladder,
                linear=not data.circular_sequences,
                min_bands=mini, max_bands=maxi,
                max_enzymes_per_digestion=data.max_enzymes
            )
        else:
            problem = SeparatingDigestionsProblem(
                 sequences=sequences, enzymes=enzymes, ladder=ladder,
                 linear=not data.circular_sequences,
                 max_enzymes_per_digestion=data.max_enzymes
            )
        self.logger(message='Selecting digestions...')
        score, selected_digestions = problem.select_digestions(
            max_digestions=data.max_digestions, search='full')
        bands_props = None if not data.show_bands_sizes else dict(
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
        
        if data.plot_cuts:
            ladder = bandwagon.custom_ladder(None, ladder.bands)
            self.logger(message="Plotting cuts maps...")
            zip_root = flametree.file_tree("@memory")
            bandwagon.plot_records_digestions(
                target=zip_root._file("Details.pdf").open("wb"),
                ladder=ladder,
                records_and_digestions=[
                    (rec, digestion)
                    for rec in records
                    for digestion in selected_digestions
                ]
            )
            pdf_data = zip_root["Details.pdf"].read('rb')
            pdf_data = data_to_html_data(pdf_data, datatype='pdf')
        else:
            pdf_data = None

        return {
            'figure_data': figure_data,
            'digestions': selected_digestions,
            'score': score,
            'pdf_data': pdf_data
        }

class SelectDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
