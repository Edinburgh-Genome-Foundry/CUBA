from .base import SerializerView
from rest_framework import status, serializers
from rest_framework.response import Response
from rest_framework.renderers import JSONRenderer
import zipsnfolders
import base64

try:
    from weasyprint import HTML
    REPORTS_AVAILABLE = True
except ImportError:
    REPORTS_AVAILABLE = False

class ReportsView(SerializerView):
    renderer_classes = (JSONRenderer,)

    class serializer_class(serializers.Serializer):
        text = serializers.CharField(label="The text for the report")
        generate_report = serializers.BooleanField(label="PDF/Zip report ?")
        generate_preview = serializers.BooleanField(label="Html preview ?")

    def post(self, request, format=None):
        """ Some description for posts"""
        data = self.serialize(request)
        if not isinstance(data, dict):
            return data
        weasywriter = HTML(string='<h1>Report</h1><p>%s</p>' % data.text)
        zip_data = zipsnfolders.write_files(target="@memory", files_contents=[
            ("README.txt", "Thank you for ordering your new pet with us !"),
            ("Report.pdf", weasywriter.write_pdf)
        ])
        response = {}
        print(base64.b64encode(zip_data)[:50])
        if data.generate_report:
            response["file"] = {
                "name": "report.zip",
                "data": 'data:application/zip;base64,' + base64.b64encode(zip_data).decode("utf-8"),
                "mimetype": "application/zip"
            }
        if data.generate_preview:
            response["preview"] = {
                "html": "<p>%s<p>" % data.text
            }
        return Response(response, status=status.HTTP_200_OK)
