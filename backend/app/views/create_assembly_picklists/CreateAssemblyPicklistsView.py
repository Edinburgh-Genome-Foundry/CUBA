"""Bla."""

from rest_framework import serializers
from plateo import AssemblyPlan
import pandas
from plateo.parsers import plate_from_content_spreadsheet
from plateo.containers.plates import Plate4ti0960
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
    source_plate = FileSerializer(allow_null=False)
    picklist = FileSerializer()
    destination_type = serializers.CharField()
    destination_size = serializers.IntegerField()
    fill_by = serializers.CharField()
    destination_plate = FileSerializer(allow_null=True)
    quantity_unit = serializers.CharField()
    part_quantity = serializers.FloatField()
    buffer_volume = serializers.FloatField()
    total_volume = serializers.FloatField()
    parts_infos = serializers.ListField(child=FileSerializer())
    dispenser_machine = serializers.CharField()
    dispenser_min_volume = serializers.FloatField()
    dispenser_max_volume = serializers.FloatField()
    dispenser_resolution = serializers.FloatField()


class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data

        # root = flametree.file_tree(".")
        picklist_filelike = file_to_filelike_object(data.picklist)
        if data.picklist.name.endswith('.csv'):
            dataframe = pandas.read_csv(picklist_filelike, index=0)
        else:
            dataframe = pandas.read_excel(picklist_filelike, index=0)
        assembly_plan = AssemblyPlan(OrderedDict([
            (row[0], [e for e in row[1:] if str(e) not in ['-', 'nan']])
            for i, row in dataframe.iterrows()
            if row[0] not in ['nan', 'Construct name']
        ]))
        for assembly, parts in assembly_plan.assemblies.items():
            assembly_plan.assemblies[assembly] = [
                part.replace(" ", "_") for part in parts
            ]

        if len(data.parts_infos):
            first_file = data.parts_infos[0]
            if first_file.name.endswith(('.csv', '.xls', '.xlsx')):
                first_file_filelike = file_to_filelike_object(first_file)
                if first_file.name.endswith('.csv'):
                    dataframe = pandas.read_csv(first_file_filelike)
                else:
                    dataframe = pandas.read_excel(first_file_filelike)
                parts_data = {
                    row.part: {'size': row['size']}
                    for i, row in dataframe.iterrows()
                }
            else:
                records = records_from_data_files(data.parts_infos)
                parts_data = {
                    rec.file_name: {'record': rec}
                    for rec in records
                }
            assembly_plan.parts_data = parts_data
            parts_without_data = assembly_plan.parts_without_data()
            if len(parts_without_data):
                return {
                    'success': False,
                    'message': 'Some parts have no provided record or data.',
                    'missing_parts': parts_without_data
                }

        if data.quantity_unit == 'fmol':
            part_mol = data.part_quantity * 1e-15
            part_g = None
        if data.quantity_unit == 'nM':
            part_mol = data.part_quantity * data.total_volume * 1e-15
            part_g = None
        if data.quantity_unit == 'fmol':
            part_mol = None
            part_g = data.part_quantity * 1e-9

        picklist_generator = AssemblyPicklistGenerator(
            part_mol=part_mol,
            part_g=part_g,
            complement_to=data.total_volume * 1e-6,
            buffer_volume=data.buffer_volume * 1e-6,
            volume_rounding=2.5e-9,
            minimal_dispense_volume=5e-9
        )
        source_filelike = file_to_filelike_object(data.source_plate)
        source_plate = plate_from_content_spreadsheet(source_filelike)
        for well in source_plate.iter_wells():
            if well.is_empty:
                continue
            quantities = well.content.quantities
            part, quantity = quantities.items()[0]
            quantities.pop(part)
            quantities[part.replace(" ", "_")] = quantity

        print (source_plate)
        parts_on_plate = set([
            well.content.components_as_string()
            for well in source_plate.iter_wells()
            ]
        ).difference({''})
        print (parts_on_plate)

        source_plate.name = "Source"
        for well in source_plate.iter_wells():
            if not well.is_empty:
                content = well.content.components_as_string()
                well.content.quantities[content] *= 1e-3

        self.logger(message="Generating Picklist...")
        destination_plate = Plate4ti0960("Mixplate")
        picklist, picklist_data = picklist_generator.make_picklist(
            assembly_plan,
            source_wells=source_plate.iter_wells(),
            destination_wells=destination_plate.iter_wells(direction='column'),
            # complement_well=source_plate.wells.O24,
            # buffer_well=source_plate.wells.P24
        )
        if picklist is None:
            return {
                'success': False,
                'message': 'Some parts in the assembly plan have no '
                           'corresponding well.',
                'picklist_data': picklist_data,
                'missing_parts': picklist_data.get('missing_parts', None)
            }
        future_plates = picklist.execute(inplace=False)

        def text(w):
            txt = human_volume(w.content.volume)
            if 'construct' in w.data:
                txt = "\n".join([w.data.construct, txt])
            return txt
        plotter = PlateTextPlotter(text)
        ax, _ = plotter.plot_plate(future_plates[destination_plate], figsize=(20, 8))

        ziproot = flametree.file_tree("@memory", replace=True)

        # MIXPLATE MAP PLOT

        ax.figure.savefig(
            ziproot._file("final_mixplate.pdf").open('wb'),
            format="pdf",
            bbox_inches="tight")
        plt.close(ax.figure)

        self.logger(message="Writing report...")

        # ASSEMBLY REPORT

        picklist_to_assembly_mix_report(
            picklist,
            ziproot._file("assembly_mix_picklist_report.pdf").open('wb'),
            data=picklist_data)
        assembly_plan.write_report(
            ziproot._file("assembly_plan_summary.pdf").open('wb'))

        # MACHINE PICKLIST

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

class CreateAssemblyPicklistsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
