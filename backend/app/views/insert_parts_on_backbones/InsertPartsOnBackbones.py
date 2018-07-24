"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data)
from ..serializers import FileSerializer
from dnacauldron import (swap_donor_vector_part, autoselect_enzyme,
                         BackboneChoice, insert_parts_on_backbones)
import flametree
from Bio import SeqIO


digestion = serializers.ListField(child=serializers.CharField())


class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    backbone = FileSerializer(allow_null=True)
    backbones = serializers.ListField(child=FileSerializer())
    inserts = serializers.ListField(child=FileSerializer())
    mode = serializers.CharField()

class worker_class(AsyncWorker):

    def work(self):
        if self.data.mode =='one_backbone':
            return self.one_backbone()
        else:
            return self.select_backbone()

    def one_backbone(self):
        self.logger(message="Reading Data...")
        data = self.data
        backbone = records_from_data_files([data.backbone])[0]
        inserts = [
            records_from_data_files([f])[0]
            for f in data.inserts
        ]
        records = inserts + [backbone]
        for record in records:
            record.linear = False  # Trick
        if data.enzyme == "Autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)
        zip_root = flametree.file_tree('@memory')
        for insert in inserts:
            record = swap_donor_vector_part(backbone, insert, data.enzyme)
            record.id = insert.id
            print ("record id", record.id)
            SeqIO.write(record, zip_root._file(record.id + '.gb'), 'genbank')

        if len(inserts) == 1:
            f = zip_root._all_files[0]
            data = f.read('rb')
            print (zip_root._all_files)
            print ("DAAAATAAAA", data)
            return {
              'file': {
                  'data': data_to_html_data(data, 'genbank'),
                  'name': f._name,
                  'mimetype': 'application/genbank'
              },
              'success': 'true',
              'summary': 'Swapping succesful !'
            }
        else:
            return {
              'file': {
                  'data': data_to_html_data(zip_root._close(), 'zip'),
                  'name': 'donor_swap_genbanks.zip',
                  'mimetype': 'application/zip'
              },
              'success': 'yeah!',
              'summary': 'none yet'
            }

    def select_backbone(self):
        self.logger(message="Reading Data...")
        data = self.data
        backbones = records_from_data_files(data.backbones)
        inserts = records_from_data_files(data.inserts)
        # inserts = [
        #     records_from_data_files([f])[0]
        #     for f in data.inserts
        # ]
        records = inserts + backbones
        for record in records:
            record.linear = False  # Trick
        if data.enzyme == "Autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)
        zip_root = flametree.file_tree('@memory')
        choices = insert_parts_on_backbones(
            inserts, backbones, process_parts_with_backbone=True)
        dataframe = BackboneChoice.list_to_infos_spreadsheet(choices)
        dataframe.to_excel(zip_root._file('summary.xls').open('wb'), index=False)
        BackboneChoice.write_final_records(
            choices, zip_root._dir("records")._path)

        return {
          'file': {
              'data': data_to_html_data(zip_root._close(), 'zip'),
              'name': 'backbone_autoselection.zip',
              'mimetype': 'application/zip'
          },
          'success': 'yeah!',
          'summary': 'none yet'
        }

class InsertPartsOnBackbonesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
