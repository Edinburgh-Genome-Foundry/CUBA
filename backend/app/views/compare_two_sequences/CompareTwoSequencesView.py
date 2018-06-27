"""Bla."""

from rest_framework import serializers
from Bio import SeqIO
import flametree
from geneblocks import DiffBlocks

from ..base import AsyncWorker, StartJobView
from ..serializers import FileSerializer
from ..tools import (records_from_data_files, data_to_html_data,
                     record_to_formated_string,
                     matplotlib_figure_to_svg_base64_data,
                     figures_to_pdf_report_data)



class serializer_class(serializers.Serializer):
    """Serializer."""
    sequence1 = FileSerializer()
    sequence2 = FileSerializer()
    figure_width = serializers.IntegerField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data


        self.logger(message='Reading the files...')
        seq_1 = records_from_data_files([data.sequence1])[0]
        seq_2 = records_from_data_files([data.sequence2])[0]

        self.logger(message='Computing the difference blocks...')
        diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)

        self.logger(message='Computing the difference blocks...')
        ax = diff_blocks.plot(figure_width=data.figure_width)
        if not hasattr(ax, 'figure'):
            ax = ax[0]
        ax.set_title("%s, with annotated diffs to %s" %
                     (seq_2.name, seq_1.name))
        figure_data = matplotlib_figure_to_svg_base64_data(
            ax.figure, bbox_inches="tight")

        diff_features = diff_blocks.diffs_as_features()
        for f in diff_features:
            f.type = "misc_feature"
        seq_2.features += diff_features
        genbank_data = record_to_formated_string(seq_2)

        return {
          'record': {
              'data': genbank_data,
              'name': 'diff.gb',
              'mimetype': 'application/genbank'
          },
          'figure_data': figure_data
        }




class CompareTwoSequencesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
