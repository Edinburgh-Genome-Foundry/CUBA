"""Bla."""

from rest_framework import serializers
import flametree
from geneblocks import CommonBlocks

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data,
                     matplotlib_figure_to_svg_base64_data,
                     figures_to_pdf_report_data, write_record,
                     autoname_genbank_file)


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

        sequences = {
          record.name: str(record.seq).upper()
          for record in records_from_data_files(data.files)
        }

        self.logger(message='Analyzing the sequence...')
        common_blocks = CommonBlocks.from_sequences(
            sequences, block_selection_method=data.block_selection,
            min_block_size=data.min_block_size
        )


        self.logger(message='Generating the report')

        root = flametree.file_tree("@memory")

        axes = common_blocks.plot_common_blocks()
        figure = axes[0].figure
        figure_data = matplotlib_figure_to_svg_base64_data(
            figure, bbox_inches="tight")
        figure.savefig(root._file('report.pdf').open('wb'), format='pdf',
                    bbox_inches='tight')
        
        blocks = {'common': common_blocks.common_blocks_records(),
                  'unique': common_blocks.unique_blocks_records()}

        for folder, blocks_list in blocks.items():
            gb_folder = root._dir('%s_blocks_genbanks' % folder)
            for rec in blocks_list:
                gb_file = gb_folder._file(autoname_genbank_file(rec))
                write_record(rec, gb_file, 'genbank')
        
        root._file('all_blocks.fa').write("\n\n".join(
            "> %s\n%s" % (rec.id, rec.seq)
            for rec in (blocks['common'] + blocks['unique'])
        ))

        seq_gb_dir = root._dir('sequences_genbanks')
        seq_records = common_blocks.sequences_with_annotated_blocks()
        for recname, rec in seq_records.items():
            rec.id = recname
            gb_file = seq_gb_dir._file(autoname_genbank_file(rec))
            write_record(rec, gb_file, 'genbank')
        zip_data = root._close()

        return {
          'zip_file': {
              'data': data_to_html_data(zip_data, 'zip'),
              'name': 'common_blocks_report.zip',
              'mimetype': 'application/zip'
          },
          'figure_data': figure_data
        }




class FindCommonBlocksView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
