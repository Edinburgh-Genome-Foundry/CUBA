"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (data_to_html_data, records_from_data_files,
                     matplotlib_figure_to_svg_base64_data, csv_to_list)
from bandwitch import LADDERS
import bandwagon
import flametree

digestion = serializers.ListField(child=serializers.CharField())
class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class serializer_class(serializers.Serializer):
    ladder = serializers.CharField()
    digestions = serializers.ListField(
        child=serializers.ListField(
            child=serializers.CharField()))
    files = serializers.ListField(child=FileSerializer())
    circular_sequences = serializers.BooleanField()
    make_cuts_position_report = serializers.BooleanField()
    show_band_sizes = serializers.BooleanField()
    use_file_names_as_ids = serializers.BooleanField()
    use_ordering_list = serializers.BooleanField()
    ordering_list = serializers.CharField(allow_blank=True)

class worker_class(AsyncWorker):
    generate_report = True
    generate_preview = True

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data
        records = records_from_data_files(data.files)
        for f, record in zip(data.files, records):
            record.linear = not data.circular_sequences
        if data.use_file_names_as_ids:
            for r in records:
                r.id = r.name = r.file_name
        if data.use_ordering_list:
            order = csv_to_list(data.ordering_list)
            records = sorted(records, key=lambda r: order.index(r.id))
        ladder = bandwagon.custom_ladder(None, LADDERS[data.ladder].bands)
        digestions = [tuple(enzs) for enzs in data.digestions]

        self.logger(message="Generating report...")
        
        pdf_data = None
        if data.make_cuts_position_report:
            zip_root = flametree.file_tree("@memory")
            bandwagon.plot_records_digestions(
                target=zip_root._file("Details.pdf").open("wb"),
                ladder=ladder,
                records=records,
                digestions=digestions
            )
            pdf_data = zip_root["Details.pdf"].read('rb')
            pdf_data = data_to_html_data(pdf_data, datatype='pdf')
        axes = bandwagon.plot_all_digestion_patterns(
            records=records,
            digestions=digestions,
            ladder=ladder,
            group_by="digestions",
            show_band_sizes=data.show_band_sizes)
        return dict(
            success=True,
            pdf_data=pdf_data,
            figure_data=matplotlib_figure_to_svg_base64_data(
                axes[0].figure, bbox_inches='tight')
        )


class PredictDigestionsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
