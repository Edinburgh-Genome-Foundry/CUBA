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
    constructs = serializers.ListField(child=FileSerializer())
    availablePrimers = serializers.ListField(child=FileSerializer())
    circularSequences = serializers.BooleanField()

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
        selector = PrimerSelector(read_range=(150, 800), tm_range=(55, 70),
                                  logger=self.logger)
        selected_primers = selector.select_primers(
            records, available_primers)


        # PLOT THE PREDICTED SEQUENCING COVERAGE FOR EACH CONSTRUCT

        selector.plot_coverage(
            records=records,
            selected_primers=selected_primers,
            pdf_path=zip_root._file('coverages_plots.pdf').open('wb'))

        # WRITE ALL PRIMERS IN A CSV FILE (PRIMERS TO ORDER ARE FIRST)
        selector.write_primers_table(
            primers=selected_primers,
            csv_path=zip_root._file('primers_list.csv')
        )
        n_new, n_available = 1, 2
        return {
            'zip_file': {
                'data': data_to_html_data(zip_root._close(), 'zip'),
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
