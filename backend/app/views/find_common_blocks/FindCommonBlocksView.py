"""Bla."""

from rest_framework import serializers
from Bio import SeqIO
import flametree
from geneblocks import BlocksFinder

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_file, zip_data_to_html_data,
                     matplotlib_figure_to_svg_base64_data,
                     figures_to_pdf_report_data)


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


        self.logger(message='Reading the files...')

        sequences = {}
        for file_ in data.files:
            records, fmt = records_from_data_file(file_)
            single_record = len(records) == 1
            for i, record in enumerate(records):
                name = record.id
                if name in [None, '', "<unknown id>"]:
                    name = file_.name + ('' if single_record else ("%04d" % i))
                sequences[name] = str(record.seq).upper()

        self.logger(message='Analyzing the sequence...')
        blocks_finder = BlocksFinder(
            sequences, block_selection=data.block_selection,
            min_block_size=data.min_block_size
        )

        blocks_finder = BlocksFinder(sequences)


        self.logger(message='Generating the report')
        root = flametree.file_tree("@memory")

        axes = blocks_finder.plot_common_blocks()
        figure = axes[0].figure
        figure_data = matplotlib_figure_to_svg_base64_data(
            figure, bbox_inches="tight")
        figure.savefig(root._file('report.pdf').open('wb'), format='pdf',
                    bbox_inches='tight')

        genbank_dir = root._dir('blocks_genbanks')
        for record in blocks_finder.common_blocks_records():
            print (record.id)
            SeqIO.write(record, genbank_dir._file(record.id + '.gb').open('w'),
                        'genbank')
        zip_data = root._close()

        return {
          'zip_file': {
              'data': zip_data_to_html_data(zip_data),
              'name': 'optimization_report.zip',
              'mimetype': 'application/zip'
          },
          'figure_data': figure_data
        }




class FindCommonBlocksView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
