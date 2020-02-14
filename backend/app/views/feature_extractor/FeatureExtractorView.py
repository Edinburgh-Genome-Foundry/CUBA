"""Bla."""

from rest_framework import serializers

from easy_dna import extract_from_input
import flametree

from ..base import AsyncWorker, StartJobView
from ..tools import (
    records_from_data_files,
    data_to_html_data,
    matplotlib_figure_to_svg_base64_data,
    figures_to_pdf_report_data,
    write_record,
    autoname_genbank_file,
)


class FileSerializer(serializers.Serializer):
    """Serializer for files."""

    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    """Serializer."""

    files = serializers.ListField(child=FileSerializer())
    min_sequence_length = serializers.IntegerField()
    direct_sense = serializers.BooleanField()


class worker_class(AsyncWorker):
    def work(self):
        data = self.data
        direct_sense = data.direct_sense
        min_sequence_length = data.min_sequence_length

        self.logger(message="Reading the files...")

        construct_list = records_from_data_files(data.files)

        self.logger(message="Generating the report")

        root = flametree.file_tree("@memory")

        extract_from_input(
            construct_list=construct_list,
            direct_sense=direct_sense,
            output_path=root,
            min_sequence_length=min_sequence_length,
        )
        zip_data = root._close()

        return {
            "zip_file": {
                "data": data_to_html_data(zip_data, "zip"),
                "name": "extracted_features.zip",
                "mimetype": "application/zip",
            }
        }


class FeatureExtractorView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
