"""Bla."""

from base64 import b64decode, b64encode
from collections import OrderedDict

from rest_framework import serializers

import flametree

from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files, data_to_html_data)
from ..serializers import FileSerializer

from primavera import PrimerSelector, Primer, load_record

class serializer_class(serializers.Serializer):
    goal = serializers.CharField()
    strand = serializers.CharField()
    constructs = serializers.ListField(child=FileSerializer())
    availablePrimers = serializers.ListField(child=FileSerializer())
    circularSequences = serializers.BooleanField()
    readRange = serializers.ListField(child=serializers.IntegerField())
    tmRange = serializers.ListField(child=serializers.IntegerField())


class worker_class(AsyncWorker):

    def work(self):

        self.logger(message='Exploring possible digestions...')

        data = self.data
        zip_root = flametree.file_tree('@memory')

        records = records_from_data_files(data.constructs)
        primers_records = records_from_data_files(data.availablePrimers)

        available_primers = [
            Primer(name=r.id, sequence=str(r.seq).upper())
            for r in primers_records
        ]

        # SELECT THE BEST PRIMERS
        #(150, 800)
        #(55, 70)
        selector = PrimerSelector(read_range=data.readRange,
                                  tm_range=data.tmRange,
                                  nucleotide_resolution=1,
                                  logger=self.logger)
        selected_primers = selector.select_primers(
            records, available_primers, strand=data.strand)

        zip_data = selector.write_multifile_report(
            records=records,
            selected_primers=selected_primers,
            target='@memory'
        )
        df = selector.write_primers_table(
            selected_primers=selected_primers,
            csv_path=zip_root._file('primers_list.csv')
        )

        n_available = df['available'].sum()
        n_new = len(df) - n_available

        return {
            'zip_file': {
                'data': data_to_html_data(zip_data, 'zip'),
                'name': 'selected_primers.zip',
                'mimetype': 'application/zip'
            },
            'success': 'yeah!',
            'summary': ('The solution found involves %d new primers and %s '
                        'available primers') % (n_new, n_available)
        }

class SelectPrimersView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
