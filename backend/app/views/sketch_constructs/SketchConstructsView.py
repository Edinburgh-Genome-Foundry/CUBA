"""Bla."""

from base64 import b64encode
import json
from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..serializers import FileSerializer
from caravagene import ConstructList


digestion = serializers.ListField(child=serializers.CharField())
class SequenceFileSerializer( FileSerializer):
    circularity = serializers.BooleanField()

class serializer_class(serializers.Serializer):
    format = serializers.CharField()
    orientation = serializers.CharField()
    fontsize = serializers.IntegerField()
    font = serializers.CharField()
    sketchesData = serializers.JSONField()

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Rendering...")
        data = self.data
        constructs_dict = data.sketchesData
        constructs_dict.update(dict(
            font=data.font,
            orientation=data.orientation,
            size=data.fontsize
        ))
        constructs = ConstructList.from_dict(constructs_dict)

        if data.format == 'PDF':
            pdf_data = constructs.to_pdf(outfile=None)
            file_data = ('data:application/pdf;base64,' +
                         b64encode(pdf_data).decode("utf-8"))
            file_mimetype = 'application/pdf'
            file_name = constructs.title + '.pdf'
        else:
            ext = format.lower()
            img_data = constructs.to_image(outfile=None, extension=ext)
            file_data = ('data:application/%s;base64,' % ext +
                         b64encode(img_data).decode("utf-8"))
            file_mimetype = 'application/' + ext
            file_name = constructs.title + '.' + ext
        return JobResult(
            file_data=file_data,
            file_mimetype=file_mimetype,
            file_name=file_name
        )

class SketchConstructsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
