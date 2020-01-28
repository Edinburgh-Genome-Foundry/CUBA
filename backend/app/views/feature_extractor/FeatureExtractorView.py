"""Bla."""

from rest_framework import serializers

import easy_dna
from easy_dna.extractor import extract_features, make_part_dict, process_report
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

        records_dict = dict()
        for input_record in construct_list:
            records = extract_features(input_record, direct_sense=direct_sense)
            key = input_record.name[0:20]  # GenBank format hard limit for name
            records_dict[key] = records

        parts_report = make_part_dict(
            records_dict, min_sequence_length=min_sequence_length
        )
        common_parts_dict = parts_report[0]
        records_dict["common_parts"] = list(common_parts_dict.values())

        for key, v in records_dict.items():
            records = records_dict[key]
            record_dir = root._dir(key)
            for j, record in enumerate(records):

                record_name_alnum = "".join(
                    x if x.isalnum() else "_" for x in record.name
                )
                record_filename = record_name_alnum + ".gb"
                record_file_path = record_dir._file(record_filename)

                try:
                    write_record(record, record_file_path, fmt="genbank")

                except Exception as err:
                    print("Error writing", record_filename, str(err))

        processed_report = process_report(parts_report[1])
        processed_report.to_csv(root._file("report.csv").open("w"))
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
