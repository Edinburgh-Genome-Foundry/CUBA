"""Bla."""

from base64 import b64encode
import json
from rest_framework import serializers
from ..base import AsyncWorker, StartJobView, JobResult
from ..serializers import FileSerializer
from caravagene import ConstructList

class serializer_class(serializers.Serializer):
    format = serializers.CharField()
    orientation = serializers.CharField()
    fontsize = serializers.IntegerField()
    font = serializers.CharField()
    width = serializers.IntegerField()
    sketchesData = serializers.JSONField()

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Rendering...")
        data = self.data
        constructs_dict = data.sketchesData
        constructs_dict.update(dict(
            font=data.font,
            orientation=data.orientation,
            size=data.fontsize,
            width=data.width
        ))
        constructs = ConstructList.from_dict(constructs_dict)
        basename = 'plot' if (constructs.title == '') else constructs.title

        if data.format == 'PDF':
            pdf_data = constructs.to_pdf(outfile=None)
            file_data = ('data:application/pdf;base64,' +
                         b64encode(pdf_data).decode("utf-8"))
            file_mimetype = 'application/pdf'
            file_name = basename + '.pdf'
        else:
            ext = data.format.lower()
            img_data = constructs.to_image(outfile=None, extension=ext)
            file_data = ('data:application/%s;base64,' % ext +
                         b64encode(img_data).decode("utf-8"))
            file_mimetype = 'application/' + ext
            file_name = basename + '.' + ext
        return JobResult(
            file_data=file_data,
            file_mimetype=file_mimetype,
            file_name=file_name
        )

class SketchConstructsView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
