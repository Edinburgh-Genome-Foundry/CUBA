"""Bla."""

from rest_framework import serializers
from saboteurs import (
    csv_to_groups_data,
    find_statistical_saboteurs,
    find_logical_saboteurs,
    statistics_report,
)
from ..base import AsyncWorker, StartJobView
from ..tools import data_to_html_data, file_to_filelike_object
from ..serializers import FileSerializer


class serializer_class(serializers.Serializer):
    """Serializer."""

    assemblies_data_file = FileSerializer()
    method = serializers.CharField()
    parts_input_type = serializers.CharField()
    constructs = serializers.ListField(
        child=FileSerializer(), allow_null=True, required=False
    )


class worker_class(AsyncWorker):
    def work(self):
        data = self.data
        self.logger(message="analyzing results...")
        asm_data = file_to_filelike_object(data.assemblies_data_file).read()
        if data.method == "statistical":
            groups_data = csv_to_groups_data(csv_string=asm_data.decode())
            analysis_results = find_statistical_saboteurs(groups_data)
            self.logger(message="Writing a report...")
            report_data = statistics_report(
                analysis_results,
                "@memory",
                replacements=[
                    ("groups", "assemblies"),
                    ("group", "assembly"),
                    ("member", "part"),
                ],
            )
            return {
                "report": {
                    "data": data_to_html_data(report_data, "pdf"),
                    "name": "saboteur_report.pdf",
                    "mimetype": "application/pdf",
                }
            }
        else:
            groups_data = csv_to_groups_data(csv_string=asm_data.decode())
            groups = {
                e['id']: e['members']
                for e in groups_data.values()
            }
            failed = [
                e['id']
                for e in groups_data.values()
                if e['attempts'] == e['failures']
            ]
            result = find_logical_saboteurs(
                groups, failed_groups=failed)
            return dict(
                saboteurs=result["saboteurs"], suspicious=result["suspicious"]
            )


class FindSaboteurPartsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
