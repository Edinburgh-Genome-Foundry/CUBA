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
    parts = serializers.ListField(child=FileSerializer())
    connectors = serializers.ListField(child=FileSerializer())
    include_fragment_plots = serializers.CharField()
    include_graph_plots = serializers.CharField()
    include_assembly_plots = serializers.BooleanField()
    use_file_names_as_ids = serializers.BooleanField()
    assembly_plan = FileSerializer(allow_null=True)
    topology = serializers.CharField()


class worker_class(AsyncWorker):
    def work(self):
        data = self.data
        logger = self.logger
        logger(message="Reading the Data...")

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

        # INITIALIZE ALL RECORDS IN A SEQUENCE REPOSITORY

        logger(message="Parsing the sequences...")

        records = records_from_data_files(
            data.parts, use_file_names_as_ids=data.use_file_names_as_ids
        )
        repository = dc.SequenceRepository()
        for record in records:
            # location-less features can cause bug when concatenating records.
            record.features = [
                f
                for f in record.features
                if f.location is not None
                and f.location.start <= f.location.end
            ]
        for r in records:
            set_record_topology(r, topology=data.topology)
            r.seq = r.seq.upper()
        repository.add_records(records, collection="parts")

        # CREATE A CONNECTORS COLLECTION IF CONNECTORS ARE PROVIDED
        if len(data.connectors):
            connector_records = records_from_data_files(data.connectors)
            for r in connector_records:
                set_record_topology(r, topology=data.topology)
                r.seq = r.seq.upper()
                if hasattr(r, "zip_file_name"):
                    collection_name = ".".join(r.zip_file_name.split(".")[:-1])
                else:
                    collection_name = ".".join(r.file_name.split(".")[:-1])
                repository.add_record(r, collection=collection_name)

        logger(message="Simulating the assembly plan...")
        filelike = file_to_filelike_object(data.assembly_plan)
        assembly_plan = dc.AssemblyPlan.from_spreadsheet(
            assembly_class="from_spreadsheet",
            path=filelike,
            name="_".join(data.assembly_plan.name.split(".")[:-1]),
            is_csv=data.assembly_plan.name.lower().endswith(".csv"),
            logger=logger,
        )
        simulation = assembly_plan.simulate(sequence_repository=repository)
        stats = simulation.compute_stats()
        n_errors = stats["errored_assemblies"]
        logger(submessage="%s error(s) found" % n_errors)
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
        logger(message="All done!")

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


class SimulateMultiMethodAssembliesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
