"""Bla."""

from matplotlib.backends.backend_pdf import PdfPages

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files, record_to_formated_string,
                     data_to_html_data)
from io import BytesIO
import flametree

class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class serializer_class(serializers.Serializer):
    format = serializers.CharField()
    inSingleFile = serializers.BooleanField()
    files = serializers.ListField(child=FileSerializer())
    use_file_names_as_ids = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):

        data = self.data
        fmt = data.format
        records = records_from_data_files(data.files)
        if data.use_file_names_as_ids:
            for r in records:
                r.id = r.name = r.file_name
        for r in records:
            r.id = r.name = r.id.replace(' ', '_')
        extension = {'fasta': '.fa', 'genbank': '.gb', 'csv': '.csv'}[fmt]
        if fmt == 'csv':
            file_content = ("name,sequence\n" + "\n".join([
                ",".join([r.id, str(r.seq)])
                for r in records
            ]))
            return JobResult(
                file_name='sequences.csv',
                file_data=file_content,
                file_mimetype='application/csv'
            )
        if (fmt == 'fasta') and (len(records) > 1) and data.inSingleFile:
            file_content = record_to_formated_string(records, fmt=fmt)
            file_content = file_content.decode().replace(' .\n', '\n').encode()
            return JobResult(
                file_name='sequences.fasta',
                file_data=data_to_html_data(bytes(file_content), 'fasta'),
                file_mimetype='application/fasta'
            )
        elif (len(records) == 1):
            record = records[0]
            file_content = record_to_formated_string(record, fmt=fmt)
            return JobResult(
                file_name=record.id + extension,
                file_data=data_to_html_data(file_content, fmt),
                file_mimetype='application/' + fmt
            )
        else:
            zip_root = flametree.file_tree('@memory')
            for record in records:
                file_content = record_to_formated_string(record, fmt=fmt)
                zip_root._file(record.id + extension).write(file_content)
            return JobResult(
                file_name='sequences.zip',
                file_data=data_to_html_data(zip_root._close(), 'zip'),
                file_mimetype='application/' + fmt
            )

class ConvertSequenceFilesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
