"""Bla."""

from collections import OrderedDict

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data,
                     file_to_filelike_object)
from ..serializers import FileSerializer

import flametree
from plateo.parsers import plate_from_platemap_spreadsheet
from bandwitch import Clone, BandsObservation, ClonesObservations


digestion = serializers.ListField(child=serializers.CharField())
class SequenceFileSerializer( FileSerializer):
    circularity = serializers.BooleanField()

class serializer_class(serializers.Serializer):
    constructsMap = FileSerializer(allow_null=True)
    clonesMap = FileSerializer(allow_null=True)
    constructsSequences = serializers.ListField(child=FileSerializer())
    goal = serializers.CharField()
    uniqueDigestion = serializers.BooleanField()
    digestion = serializers.ListField(child=serializers.CharField())
    digestionsMap = FileSerializer(allow_null=True)
    tolerance = serializers.FloatField()
    bandsRange = serializers.ListField(child=serializers.IntegerField())
    fragmentAnalysisArchive = FileSerializer(allow_null=True)
    includeDigestionPlots = serializers.BooleanField()

def file_type(f):
    return 'csv' if f.name.lower().endswith('csv') else 'excel'

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading the files...")
        data = self.data

        # PARSE ALL FILES

        constructs_records = records_from_data_files(data.constructsSequences)
        constructs_records = {r.id: r for r in constructs_records}
        print (list(constructs_records.keys()))
        constructs_records = OrderedDict(sorted(constructs_records.items()))

        constructs_plate = plate_from_platemap_spreadsheet(
            file_to_filelike_object(data.constructsMap),
            file_type=file_type(data.constructsMap),
            data_field='construct', headers=True)
        constructs_map = OrderedDict([
            (well.name, well.data.construct)
            for well in constructs_plate.iter_wells(direction='row')
            if 'construct' in well.data
            and str(well.data.construct) != 'nan'
        ])

        if data.uniqueDigestion:
            digestion = tuple(data.digestion)
            digestions_map = OrderedDict([
                (wellname, digestion)
                for wellname in constructs_map
            ])
        else:
            digestions_plate = plate_from_platemap_spreadsheet(
                file_to_filelike_object(data.digestionsMap),
                file_type=file_type(data.digestionsMap),
                data_field='digestion', headers=True)
            digestions_map = OrderedDict([
                (well.name, well.data.digestion.split(', '))
                for well in digestions_plate.iter_wells(direction='row')
                if 'digestion' in well.data
                and str(well.data.digestion) != 'nan'
            ])

        archive = file_to_filelike_object(data.fragmentAnalysisArchive)

        # ANALYZE ALL FILES AND VALIDATE BANDS

        self.logger(message="Analyzing the data...")

        observations = BandsObservation.from_aati_fa_archive(archive)
        clones = Clone.from_bands_observations(observations, constructs_map,
                                               digestions_map)
        clones_observations = ClonesObservations(clones, constructs_records)
        validations = clones_observations.validate_all_clones(
            min_band_cutoff=data.bandsRange[0],
            max_band_cutoff=data.bandsRange[1],
            relative_tolerance=data.tolerance
        )

        # CREATE A ZIP WITH VALIDATION REPORTS

        zip_root = flametree.file_tree('@memory')
        self.logger(message="Generating the validation report...")
        zip_root._file('validations.pdf').write(
            clones_observations.plot_all_validations_patterns(validations)
        )
        if data.includeDigestionPlots:
            self.logger(message="Generating the digestion sites report...")
            zip_root._file('digestions.pdf').write(
                clones_observations.plot_all_constructs_digestions()
            )

        self.logger(message="Generating the success plate map...")
        ax = clones_observations.plot_validations_plate_map(validations)
        ax.figure.savefig(zip_root._file('success_map.pdf').open('wb'),
                          format='pdf', bbox_inches='tight')

        self.logger(message="All done !")

        return {
          'zip_file': {
              'data': data_to_html_data(zip_root._close(), 'zip'),
              'name': 'validation_report.zip',
              'mimetype': 'application/zip'
          },
          'success': 'yeah!',
          'summary': 'none yet'
        }

class AnalyzeDigestsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
