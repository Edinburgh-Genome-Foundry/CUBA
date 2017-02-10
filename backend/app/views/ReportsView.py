from .base import SerializerView
from rest_framework import status, serializers
from rest_framework.response import Response
from rest_framework.renderers import JSONRenderer
import flametree
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
        print ("rolollolol", data.text)
        weasywriter = HTML(string=u"<h1>Report</h1><p>%s</p>" % data.text)
        root = flametree.file_tree("@memory")
        weasywriter.write_pdf(root._file("Report.pdf"))
        root._file("README.txt").write(
            data.text + " Thanks for ordering your pet with us!")
        zip_data = root._close()

        response = {}
        print(base64.b64encode(zip_data)[:50])
        if data.generate_report:
            response["file"] = {
                "name": "report.zip",
                "data": ('data:application/zip;base64,' +
                         base64.b64encode(zip_data).decode("utf-8")),
                "mimetype": "application/zip"
            }
        if data.generate_preview:
            response["preview"] = {
                "html": "<p>%s<p>" % data.text
            }
        return Response(response, status=status.HTTP_200_OK)
