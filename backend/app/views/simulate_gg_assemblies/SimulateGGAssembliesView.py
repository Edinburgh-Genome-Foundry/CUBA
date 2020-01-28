"""Bla."""

from base64 import b64decode, b64encode

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import (
    records_from_data_files,
    file_to_filelike_object,
    data_to_html_data,
    set_record_topology,
)
from ..serializers import FileSerializer
from dnacauldron import (
    full_assembly_report,
    autoselect_enzyme,
    full_assembly_plan_report,
)
from plateo import AssemblyPlan
import pandas

digestion = serializers.ListField(child=serializers.CharField())
# class SequenceFileSerializer(FileSerializer):
#     circularity = serializers.BooleanField()


class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    parts = serializers.ListField(child=FileSerializer())
    connectors = serializers.ListField(child=FileSerializer())
    include_fragments = serializers.CharField()
    use_assembly_plan = serializers.BooleanField()
    single_assemblies = serializers.BooleanField()
    use_file_names_as_ids = serializers.BooleanField()
    assembly_plan = FileSerializer(allow_null=True)
    show_overhangs = serializers.BooleanField()
    backbone_first = serializers.BooleanField()
    backbone_name = serializers.CharField(allow_blank=True)
    no_skipped_parts = serializers.BooleanField()
    topology = serializers.CharField()


class worker_class(AsyncWorker):
    def work(self):
        self.logger(message="Reading Data...")
        data = self.data

        records = records_from_data_files(data.parts)
        for record in records:
            # location-less features can cause bug when concatenating records.
            record.features = [
                f
                for f in record.features
                if f.location is not None
                and f.location.start <= f.location.end
            ]
        if data.use_file_names_as_ids:
            for r in records:
                r.id = r.name = r.file_name
                if data.backbone_first and r.id == data.backbone_name:
                    r.is_backbone = True
        connector_records = records_from_data_files(data.connectors)
        for r in records + connector_records:
            set_record_topology(r, topology=data.topology)
            r.seq = r.seq.upper()
        if data.enzyme == "autoselect":
            possible_enzymes = ["BsaI", "BsmBI", "BbsI", "SapI"]
            data.enzyme = autoselect_enzyme(records, enzymes=possible_enzymes)

        self.logger(message="Generating a report, be patient.")

        if data.use_assembly_plan:
            filelike = file_to_filelike_object(data.assembly_plan)
            if data.assembly_plan.name.lower().endswith(".csv"):
                content = filelike.read().decode()
                dataframe = pandas.DataFrame(
                    [
                        [e.strip() for e in line.split(",") if len(e.strip())]
                        for line in content.split("\n")
                        if len(line)
                    ]
                )
            else:

                dataframe = pandas.read_excel(filelike, header=None)
            assembly_plan = AssemblyPlan.from_spreadsheet(dataframe=dataframe)
            assembly_plan.parts_data = {r.id: {"record": r} for r in records}
            parts_without_data = assembly_plan.parts_without_data()
            if len(parts_without_data):
                return {"success": False, "unknown_parts": parts_without_data}
            data.include_fragments = {
                'yes': True,
                'no': False,
                'on_failure': 'on_failure'
            }[data.include_fragments]
            errors, zip_data = full_assembly_plan_report(
                assembly_plan.assemblies_with_records(),
                target="@memory",
                enzyme=data.enzyme,
                assert_single_assemblies=data.single_assemblies,
                logger=self.logger,
                connector_records=connector_records,
                fail_silently=False,
                include_fragments_plots=data.include_fragments,
                include_parts_plots=data.include_fragments,
                show_overhangs_in_genbank=data.show_overhangs,
                no_skipped_parts=data.no_skipped_parts,
            )
            infos = dict(errors=errors)

        else:

            nconstructs, zip_data = full_assembly_report(
                records,
                connector_records=connector_records,
                target="@memory",
                enzyme=data.enzyme,
                max_assemblies=40,
                fragments_filters="auto",
                assemblies_prefix="assembly",
                include_fragments_plots=data.include_fragments,
                include_parts_plots=data.include_fragments,
                show_overhangs_in_genbank=data.show_overhangs,
            )
            infos = dict(nconstructs=nconstructs)

        return {
            "file": {
                "data": data_to_html_data(zip_data, "zip"),
                "name": "predicted_assemblies.zip",
                "mimetype": "application/zip",
            },
            "success": True,
            "infos": infos,
        }


class SimulateGGAssembliesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class

import scipy
import networkx
print (scipy.__version__, networkx.__version__)