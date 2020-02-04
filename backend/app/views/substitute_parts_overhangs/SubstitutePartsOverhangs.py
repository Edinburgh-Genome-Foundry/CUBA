"""Bla."""

from rest_framework import serializers
from dnacauldron.utils import substitute_overhangs
from geneblocks import DiffBlocks
import flametree
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages

from ..base import AsyncWorker, StartJobView
from ..tools import (
    records_from_data_files,
    data_to_html_data,
    write_record,
    record_to_formated_string,
)
from ..serializers import FileSerializer


digestion = serializers.ListField(child=serializers.CharField())


class serializer_class(serializers.Serializer):
    parts = serializers.ListField(child=FileSerializer())
    use_file_names_as_ids = serializers.BooleanField()
    return_linear_parts = serializers.BooleanField()
    substitutions = serializers.CharField()
    enzyme = serializers.CharField()


class worker_class(AsyncWorker):
    def work(self):
        data = self.data
        records = records_from_data_files(data.parts)
        if data.use_file_names_as_ids:
            for r in records:
                r.id = r.name = r.file_name
        substitutions = {
            l.split("=>")[0].strip(): l.split("=>")[1].strip()
            for l in data.substitutions.split("\n")
            if "=>" in l
        }

        new_records = []
        for record in records:
            enzyme = data.enzyme.replace("Autoselect", "auto")
            record.annotations['topology'] = 'circular'
            new_record = substitute_overhangs(
                record,
                substitutions,
                enzyme=enzyme,
                return_linear_parts=data.return_linear_parts,
            )
            new_record.file_name = record.file_name
            new_record.id = record.id
            new_records.append(new_record)

        if data.return_linear_parts:
            data = "\n\n".join(
                [">%s\n%s" % (r.id, str(r.seq)) for r in new_records]
            )
            filedata = {
                "data": data,
                "name": "modified_records.fa",
                "mimetype": "application/biosequence.fasta",
            }
        elif len(new_records) == 1:
            data = record_to_formated_string(new_records[0])
            filedata = {
                "data": data_to_html_data(data, "txt"),
                "name": "%s_modified.gb" % new_record.file_name,
                "mimetype": "application/biosequence.genbank",
            }
        else:
            zip_root = flametree.file_tree("@memory")
            for rec in new_records:
                f = zip_root._file("%s_modified.gb" % rec.file_name)
                write_record(rec, f)
            filedata = {
                "data": data_to_html_data(zip_root._close(), "zip"),
                "name": "records_with_substitutions.zip",
                "mimetype": "application/zip",
            }
        return {
            "file": filedata,
            # "pdf_report": {
            #     "data": data_to_html_data(pdf_io.getvalue(), "pdf"),
            #     "name": "sequences_edits.pdf",
            #     "mimetype": "application/pdf",
            # },
            "success": True,
        }


class SubstitutePartsOverhangsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
