"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_file, zip_data_to_html_data,
                     matplotlib_figure_to_svg_base64_data)
from geneblocks import BlocksFinder

class FileSerializer(serializers.Serializer):
    """Serializer for files."""

    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    """Serializer."""
    files = serializers.ListField(child=FileSerializer())
    block_selection = serializers.CharField()
    min_block_size = serializers.IntegerField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data


        self.set_progress_message('Reading the files...')

        sequences = {}
        for file_ in data.files:
            records, fmt = records_from_data_file(file_)
            single_record = len(records) == 1
            for i, record in enumerate(records):
                name = record.id
                if name in [None, '', "<unknown id>"]:
                    name = file_.name + ('' if single_record else ("%04d" % i))
                sequences[name] = str(record.seq).upper()

        self.set_progress_message('Analyzing the sequence...')
        blocks_finder = BlocksFinder(
            sequences, block_selection=data.block_selection,
            min_block_size=data.min_block_size
        )

        self.set_progress_message('Generating the report')
        axes = blocks_finder.plot_common_blocks()
        figure = axes[0].figure

        figure_data = matplotlib_figure_to_svg_base64_data(
            figure, bbox_inches="tight")

        return {
          'pdf_report': None,
          'figures_data': figure_data
        }




class FindCommonBlocksView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
