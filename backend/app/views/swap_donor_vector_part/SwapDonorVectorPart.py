"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files, data_to_html_data)
from ..serializers import FileSerializer
from dnacauldron import swap_donor_vector_part, autoselect_enzyme
import flametree
from Bio import SeqIO


digestion = serializers.ListField(child=serializers.CharField())


class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    donor_vector = FileSerializer()
    inserts = serializers.ListField(child=FileSerializer())

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data
        donor_vector = records_from_data_files([data.donor_vector])[0]

        inserts = [
            records_from_data_files([f])[0]
            for f in data.inserts
        ]
        records = inserts + [donor_vector]
        for record in records:
            record.linear = False # Trick
        if data.enzyme == "Autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)
        zip_root = flametree.file_tree('@memory')
        for insert in inserts:
            record = swap_donor_vector_part(donor_vector, insert, data.enzyme)
            record.id = insert.id
            print ("record id", record.id)
            SeqIO.write(record, zip_root._file(record.id + '.gb'), 'genbank')

        if len(inserts) == 1:
            f = zip_root._all_files[0]
            data = f.read()
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

class SwapDonorVectorPartView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
