"""Bla."""

from rest_framework import serializers
from Bio import SeqIO
import flametree
from sequenticon import sequenticon_batch, sequenticon_batch_pdf

from ..base import AsyncWorker, StartJobView
from ..tools import (records_from_data_files, data_to_html_data)


class FileSerializer(serializers.Serializer):
    """Serializer for files."""

    name = serializers.CharField()
    content = serializers.CharField()


class serializer_class(serializers.Serializer):
    """Serializer."""
    files = serializers.ListField(child=FileSerializer())
    output = serializers.CharField()
    title = serializers.CharField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        records = records_from_data_files(data.files)

        if data.output == "pdf":
            pdf_data = sequenticon_batch_pdf(records, title=data.title)
            # print(pdf_data)
            return {
              'pdf_file': {
                  'data': data_to_html_data(pdf_data, 'pdf'),
                  'name': '%s.pdf' % data.title,
                  'mimetype': 'application/pdf'
              }
            }

        else:
            images_data = sequenticon_batch(records, size=120,
                                            output_format="html_image")
            return {'images': images_data}




class RenderSequenticonsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
