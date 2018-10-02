"""Bla."""

from rest_framework import serializers
from plateo import AssemblyPlan
import pandas
from plateo.parsers import plate_from_content_spreadsheet
from plateo.containers.plates import Plate4ti0960, Plate384, get_plate_class
from plateo.exporters import (AssemblyPicklistGenerator,
                              picklist_to_assembly_mix_report)
from plateo.exporters import (picklist_to_labcyte_echo_picklist_file,
                              picklist_to_tecan_evo_picklist_file,
                              PlateTextPlotter)
from plateo.tools import human_volume
import flametree
from collections import OrderedDict
import matplotlib.pyplot as plt

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, file_to_filelike_object,
                     data_to_html_data)
from ..serializers import FileSerializer


class serializer_class(serializers.Serializer):
    rearraying_type = serializers.CharField()
    sample_source_by = serializers.CharField()
    fill_destination_by = serializers.CharField()
    destination_plate = FileSerializer(allow_null=True)
    source_name = serializers.CharField()
    destination_name = serializers.CharField()
    source_plate = FileSerializer(allow_null=True)
    destination_type = serializers.CharField()
    destination_size: serializers.IntegerField()
    dispenser_machine = serializers.CharField()
    dispenser_max_volume = serializers.FloatField()
    dispenser_min_volume = serializers.FloatField()
    dispenser_resolution = serializers.FloatField()
    dispenser_dead_volume = serializers.FloatField()


class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data

        source_filelike = file_to_filelike_object(data.source_plate)
        source = plate_from_content_spreadsheet(source_filelike)
        source.name = data.source_name

        if ((data.destination_plate is not None) and
            ((data.rearraying_type == 'map') or
             (data.destination_type == 'existing'))):
            destination = plate_from_content_spreadsheet(dest_filelike)
            destination.name = destination_name
        else:
            destination = get_plate_class(data.destination_size)()
            destination.name = destination_name

        if rearraying_type == 'map':
            # for well in destination.iter_wells():
                # well.content.volume *=  1e-6
            picklist = PickList()
            for well in source.iter_wells():
                if well.is_empty:
                    continue
                part = (well.content.components_as_string())
                destination_well = destination.find_unique_well(
                    condition=lambda w: w.content.components_as_string() == part)
                picklist.add_transfer(well, destination_well, destination_well.volume)
                destination_well.empty_completely()
            picklist.execute()
            picklist_to_tecan_evo_picklist_file(picklist, "rearray_2018-10-02.gwl")
            plate_to_content_spreadsheet(destination, "destination_after_picklist.xlsx")
            


        else:
            pass

        future_plates = picklist.execute(inplace=False)

        def text(w):
            txt = human_volume(w.content.volume)
            if 'construct' in w.data:
                txt = "\n".join([w.data.construct, txt])
            return txt
        plotter = PlateTextPlotter(text)
        ax, _ = plotter.plot_plate(future_plates[destination_plate], figsize=(20, 8))

        ziproot = flametree.file_tree("@memory", replace=True)
        ax.figure.savefig(
            ziproot._file("final_mixplate.pdf").open('wb'),
            format="pdf",
            bbox_inches="tight")
        plt.close(ax.figure)
        picklist_to_assembly_mix_report(
            picklist,
            ziproot._file("assembly_mix_picklist_report.pdf").open('wb'),
            data=picklist_data)
        assembly_plan.write_report(
            ziproot._file("assembly_plan_summary.pdf").open('wb'))
        if data.dispenser_machine == 'labcyte_echo':
            picklist_to_labcyte_echo_picklist_file(
                picklist, ziproot._file("ECHO_picklist.csv").open('w'))
        else:
            picklist_to_tecan_evo_picklist_file(
                picklist, ziproot._file("EVO_picklist.gwl").open('w'))
        zip_data = ziproot._close()

        return {
             'file': {
                 'data': data_to_html_data(zip_data, 'zip'),
                 'name': 'assemblies.zip',
                 'mimetype': 'application/zip'
             },
             'success': True
        }

class RearrayPlatesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
