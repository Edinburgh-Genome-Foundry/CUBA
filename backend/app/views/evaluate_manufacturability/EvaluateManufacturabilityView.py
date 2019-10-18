"""Bla."""

from base64 import b64decode, b64encode
from matplotlib.backends.backend_pdf import PdfPages

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (
    records_from_data_file,
    data_to_html_data,
    records_from_data_files,
    matplotlib_figure_to_svg_base64_data,
    write_record,
)
from io import BytesIO, StringIO

import dnachisel as dc
import dnachisel.reports.constraints_reports as cr
import flametree

import matplotlib.pyplot as plt


class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    files = serializers.ListField(child=FileSerializer())
    show_features = serializers.BooleanField()


class worker_class(AsyncWorker):
    def work(self):

        data = self.data
        figures = []

        self.logger(message="Generating report...")
        records = records_from_data_files(data.files)
        constraints = [
            dc.AvoidPattern("BsaI_site"),
            dc.AvoidPattern("BsmBI_site"),
            dc.AvoidPattern("BbsI_site"),
            dc.AvoidPattern("8x1mer"),
            dc.AvoidPattern("5x3mer"),
            dc.AvoidPattern("9x2mer"),
            dc.AvoidHairpins(stem_size=20, hairpin_window=200),
            dc.EnforceGCContent(mini=0.3, maxi=0.7, window=100),
        ]

        dataframe = cr.constraints_breaches_dataframe(constraints, records)
        spreadsheet_io = BytesIO()
        dataframe.to_excel(spreadsheet_io)
        records = cr.records_from_breaches_dataframe(dataframe, records)
        pdf_io = BytesIO()
        cr.breaches_records_to_pdf(records, pdf_io, logger=self.logger)
        zipped_records = flametree.file_tree("@memory")
        for record in records:
            write_record(record, zipped_records._file("%s.gb" % record.id))

        return {
            "pdf_report": {
                "data": data_to_html_data(
                    pdf_io.getvalue(),
                    "pdf",
                    filename="manufacturability_report.pdf",
                ),
                "name": "manufacturability_report.pdf",
                "mimetype": "application/pdf",
            },
            "records": {
                "data": data_to_html_data(
                    zipped_records._close(),
                    "zip",
                    filename="manufacturability_annotated_records.zip",
                ),
                "name": "manufacturability_annotated_records.zip",
                "mimetype": "application/zip",
            },
            "spreadsheet": {
                "data": data_to_html_data(
                    spreadsheet_io.getvalue(),
                    "xlsx",
                    filename="manufacturability_report.xlsx",
                ),
                "name": "manufacturability_report.xlsx",
                "mimetype": "vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            },
        }


class EvaluateManufacturabilityView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
