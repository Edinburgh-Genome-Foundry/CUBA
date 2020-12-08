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
import dnacauldron as dc
import pandas


class serializer_class(serializers.Serializer):
    enzyme = serializers.CharField()
    parts = serializers.ListField(child=FileSerializer())
    connectors = serializers.ListField(child=FileSerializer())
    include_fragment_plots = serializers.CharField()
    include_graph_plots = serializers.CharField()
    include_assembly_plots = serializers.BooleanField()
    use_assembly_plan = serializers.BooleanField()
    single_assemblies = serializers.BooleanField()
    use_file_names_as_ids = serializers.BooleanField()
    assembly_plan = FileSerializer(allow_null=True)
    backbone_first = serializers.BooleanField()
    backbone_name = serializers.CharField(allow_blank=True)
    no_skipped_parts = serializers.BooleanField()
    topology = serializers.CharField()


class worker_class(AsyncWorker):
    def work(self):
        self.logger(message="Reading the Data...")
        data = self.data

        # CHANGE THE VALUE OF THE BOOLEANS

        def yes_no(value):
            return {"yes": True, "no": False}.get(value, value)

        include_fragment_plots = yes_no(data.include_fragment_plots)
        include_graph_plots = yes_no(data.include_graph_plots)
        include_assembly_plots = yes_no(data.include_assembly_plots)
        report_writer = dc.AssemblyReportWriter(
            include_mix_graphs=include_graph_plots,
            include_assembly_plots=include_assembly_plots,
            include_fragment_plots=include_fragment_plots,
        )

        # INITIALIZE ALL RECORDS IN A SEQUENCE REPOSITORY

        records = records_from_data_files(
            data.parts, use_file_names_as_ids=data.use_file_names_as_ids
        )
        repository = dc.SequenceRepository()
        for record in records:
            # location-less features can cause bug when concatenating records.
            record.features = [
                f
                for f in record.features
                if f.location is not None and f.location.start <= f.location.end
            ]
        for r in records:
            if data.backbone_first and r.id == data.backbone_name:
                r.is_backbone = True
        repository.add_records(records, collection="parts")

        # CREATE A CONNECTORS COLLECTION IF CONNECTORS ARE PROVIDED

        connectors_collection = None
        if len(data.connectors):
            connector_records = records_from_data_files(data.connectors)
            for r in records + connector_records:
                set_record_topology(r, topology=data.topology)
                r.seq = r.seq.upper()
            repository.add_records(connector_records, collection="connectors")
            connectors_collection = "connectors"

        # SIMULATE!

        self.logger(message="Simulating the assembly...")

        if not data.use_assembly_plan:

            # SCENARIO: SINGLE ASSEMBLY

            parts = [r.id for r in records]
            assembly = dc.Type2sRestrictionAssembly(
                name="simulated_assembly",
                parts=parts,
                enzyme=data.enzyme,
                expected_constructs="any_number",
                connectors_collection=connectors_collection,
            )
            simulation = assembly.simulate(sequence_repository=repository)
            n = len(simulation.construct_records)
            self.logger(message="Done (%d constructs found), writing report..." % n)
            report_zip_data = simulation.write_report(
                target="@memory", report_writer=report_writer
            )
            return {
                "file": {
                    "data": data_to_html_data(report_zip_data, "zip"),
                    "name": "assembly_simulation.zip",
                    "mimetype": "application/zip",
                },
                "errors": [str(e) for e in simulation.errors],
                "n_constructs": len(simulation.construct_records),
                "success": True,
            }

        else:

            # SCENARIO: FULL ASSEMBLY PLAN

            filelike = file_to_filelike_object(data.assembly_plan)
            assembly_plan = dc.AssemblyPlan.from_spreadsheet(
                assembly_class=dc.Type2sRestrictionAssembly,
                path=filelike,
                connectors_collection=connectors_collection,
                expect_no_unused_parts=data.no_skipped_parts,
                expected_constructs=1 if data.single_assemblies else "any_number",
                name="_".join(data.assembly_plan.name.split(".")[:-1]),
                is_csv=data.assembly_plan.name.lower().endswith(".csv"),
                logger=self.logger,
            )

            simulation = assembly_plan.simulate(sequence_repository=repository)
            stats = simulation.compute_stats()
            n_errors = stats["errored_assemblies"]
            self.logger(message="Done (%d errors), writing report..." % n_errors)
            report_zip_data = simulation.write_report(
                target="@memory",
                assembly_report_writer=report_writer,
                logger=self.logger,
            )
            errors = [
                error
                for assembly_simulation in simulation.assembly_simulations
                for error in assembly_simulation.errors
            ]

            return {
                "file": {
                    "data": data_to_html_data(report_zip_data, "zip"),
                    "name": "%s.zip" % assembly_plan.name,
                    "mimetype": "application/zip",
                },
                "assembly_stats": stats,
                "errors": [str(e) for e in errors],
                "success": True,
            }
        self.logger(message="Simulating the assembly...")


class SimulateGGAssembliesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
