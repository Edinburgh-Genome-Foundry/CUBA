"""Bla."""

from base64 import b64encode
from matplotlib.backends.backend_pdf import PdfPages

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView

from ..tools import (records_from_data_files,
                     matplotlib_figure_to_svg_base64_data)
from ..serializers import FileSerializer
from dna_features_viewer import (BiopythonTranslator, GraphicRecord,
                                 CircularGraphicRecord)
from io import BytesIO


digestion = serializers.ListField(child=serializers.CharField())

class FileSerializer(serializers.Serializer):
    name = serializers.CharField()
    content = serializers.CharField()

class StyleSerializer(serializers.Serializer):
    selector = serializers.CharField()
    feature_type = serializers.CharField()
    feature_text = serializers.CharField(allow_blank=True)
    color = serializers.CharField()
    keep_or_discard = serializers.CharField()
    display_label = serializers.BooleanField()
    thickness = serializers.IntegerField()


class serializer_class(serializers.Serializer):
    files = serializers.ListField(child=FileSerializer())
    display = serializers.CharField()
    default_color = serializers.CharField()
    default_display_label = serializers.BooleanField()
    default_thickness = serializers.IntegerField()
    plot_width = serializers.IntegerField()
    plot_ruler = serializers.BooleanField()
    inline_labels = serializers.BooleanField()
    plot_full_sequence = serializers.BooleanField()
    plot_from_position = serializers.IntegerField()
    plot_to_position = serializers.IntegerField()
    custom_styles = serializers.ListField(child=StyleSerializer())
    must_contain = serializers.CharField(allow_blank=True)
    must_not_contain = serializers.CharField(allow_blank=True)
    keep_or_discard = serializers.CharField()
    keep_or_discard_types = serializers.ListField(child=serializers.CharField())
    plot_nucleotides = serializers.BooleanField()
    plot_translation = serializers.BooleanField()
    translation_start = serializers.IntegerField()
    translation_end = serializers.IntegerField()
    pdf_report = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):
        self.logger(message="Reading Data...")
        data = self.data

        must_contain = [
            s.strip() for s in data.must_contain.split(',')
            if s.strip() != ''
        ]
        must_not_contain = [
            s.strip()
            for s in data.must_not_contain.split(',')
            if s.strip() != ''
        ]
        filter_feature_types = [f.lower() for f in data.keep_or_discard_types]

        def feature_text(f):
            return ", ".join([str(v) for v in f.qualifiers.values()])

        def feature_filter(f):
            ftype = f.type.lower()
            keep = data.keep_or_discard == 'keep'
            if filter_feature_types != []:
                in_types = ftype in filter_feature_types
                if (keep and not in_types) or (in_types and not keep):
                    return False
            text = feature_text(f)
            if len(must_contain) and not any([c in text for c in must_contain]):
                return False
            if len(must_not_contain) and any([c in text for c in must_not_contain]):
                return False
            return True

        def features_properties(f):
            properties = {
                'color': data.default_color,
                'linewidth': data.default_thickness
            }
            if not data.default_display_label:
                properties['label'] = None
            ftype = f.type.lower()
            for fl in data.custom_styles:
                keep = fl.keep_or_discard == 'keep'
                if (fl.selector == 'text'):
                    has_term = fl.feature_text in feature_text(f)
                    if (keep and has_term) or ((not keep) and (not has_term)):
                        properties['color'] = fl.color
                        properties['linewidth'] = fl.thickness
                        if fl.display_label:
                            properties.pop('label', '')
                if (fl.selector == 'type'):
                    is_type = (ftype == fl.feature_type.lower())
                    if (keep and is_type) or ((not keep) and (not is_type)):
                        properties['color'] = fl.color
                        properties['linewidth'] = fl.thickness
                        if fl.display_label:
                            properties.pop('label', '')
            return properties

        display_class = {
            'linear': GraphicRecord,
            'circular': CircularGraphicRecord
        }[data.display]

        translator = BiopythonTranslator(
            features_filters=(feature_filter,),
            features_properties=features_properties
        )
        records = records_from_data_files(data.files)
        figures = []
        for rec in self.logger.iter_bar(record=records):
            gr = translator.translate_record(rec, grecord_class=display_class)
            if not data.plot_full_sequence:
                gr = gr.crop((data.plot_from_position, data.plot_to_position))
            ax, _ = gr.plot(figure_width=data.plot_width,
                            with_ruler=data.plot_ruler,
                            annotate_inline=data.inline_labels)
            if data.plot_nucleotides:
                gr.plot_sequence(ax)
            # if data.plot_translation:
            #     gr.plot_translation
            figure = ax.figure
            figure.suptitle(rec.id)
            figures.append(figure)

        if data.pdf_report:
            pdf_io = BytesIO()

            with PdfPages(pdf_io) as pdf:
                for fig in figures:
                    pdf.savefig(fig, bbox_inches="tight")

            pdf_data = ('data:application/pdf;base64,' +
                        b64encode(pdf_io.getvalue()).decode("utf-8"))
            figures_data = {
                'data': pdf_data,
                'name': 'sequence_feature_plots.pdf',
                'mimetype': 'application/pdf'
            }
        else:
            figures_data = []
            for _file, fig in zip(data.files, figures):
                figdata = matplotlib_figure_to_svg_base64_data(
                    fig, bbox_inches="tight")
                figures_data.append({'img_data': figdata,
                                     'filename': _file.name})

        return {
          'pdf_report': None if not data.pdf_report else figures_data,
          'figures_data': None if data.pdf_report else figures_data
        }

class PlotSequenceFeaturesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
