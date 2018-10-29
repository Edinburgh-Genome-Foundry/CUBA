"""Bla."""

from rest_framework import serializers
import flametree
from geneblocks import CommonBlocks

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data,
                     matplotlib_figure_to_svg_base64_data,
                     figures_to_pdf_report_data, write_record,
                     autoname_genbank_file, write_record)


class FileSerializer(serializers.Serializer):
    """Serializer for files."""

    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    """Serializer."""
    files = serializers.ListField(child=FileSerializer())
    use_filenames_as_ids = serializers.BooleanField()
    min_block_size = serializers.IntegerField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data


        self.logger(message='Reading the files...')
        records = records_from_data_files(data.files)
        print ("loooool")
        for r in records:
            if data.use_filenames_as_ids:
                r.id = r.name = r.file_name
        

        self.logger(message='Finding transfers...')
        records_dict = {r.id: r for r in records}
        blocks = CommonBlocks(records, min_block_size=data.min_block_size)
        new_records = blocks.copy_features_between_common_blocks(inplace=False)
        added_features = [
            (name, len(new_record.features) - len(records_dict[name].features))
            for name, new_record in new_records.items()
        ]
        ziproot = flametree.file_tree('@memory')
        for name, new_record in new_records.items():
            write_record(new_record, ziproot._file(name + '.gb'))
        zip_data = ziproot._close()

        return {
          'zip_file': {
              'data': data_to_html_data(zip_data, 'zip'),
              'name': 'optimization_report.zip',
              'mimetype': 'application/zip'
          },
          'added_features': added_features
        }




class TransferFeaturesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
