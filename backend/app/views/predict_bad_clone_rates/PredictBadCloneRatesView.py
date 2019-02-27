"""Bla."""

from base64 import b64decode, b64encode
from collections import OrderedDict
from io import BytesIO
import textwrap

from rest_framework import serializers
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from plateo import AssemblyPlan
import kappagate
import pandas

from ..base import AsyncWorker, StartJobView, JobResult
from ..tools import (records_from_data_files, file_to_filelike_object,
                     data_to_html_data, csv_to_list)
from ..serializers import FileSerializer

class serializer_class(serializers.Serializer):
    input_type = serializers.CharField()
    constructs_files = serializers.ListField(child=FileSerializer())
    parts_files = serializers.ListField(child=FileSerializer())
    overhangs_list = serializers.CharField(allow_blank=True)
    use_assembly_plan = serializers.BooleanField()
    assembly_plan = FileSerializer()
    rate_limit_in_plots = serializers.IntegerField()
    use_file_names_as_ids = serializers.BooleanField()
    assembly_plan = FileSerializer(allow_null=True)
    dataset_tm = serializers.CharField()
    dataset_duration = serializers.CharField()
    dataset_corrective_factor = serializers.FloatField()
    backbone_string = serializers.CharField(allow_blank=True)

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        log = self.logger
        log(message="Reading the data...")
        if data.input_type == 'parts':

            records = records_from_data_files(data.parts_files)
            for r in records:
                if data.use_file_names_as_ids:
                    r.id = r.name = r.file_name
                if not hasattr(r, 'linear'):
                    r.linear = False
                r.seq = r.seq.upper()

            

            if data.use_assembly_plan:
                filelike = file_to_filelike_object(data.assembly_plan)
                if data.assembly_plan.name.lower().endswith('.csv'):
                    content = filelike.read().decode()
                    dataframe = pandas.DataFrame([
                        [e.strip() for e in line.split(',') if len(e.strip())]
                        for line in content.split('\n')
                        if len(line)
                    ])
                else:
                    dataframe = pandas.read_excel(filelike, header=None)
                assembly_plan = AssemblyPlan.from_spreadsheet(dataframe=dataframe)
                slots_dict = {
                    name: kappagate.parts_records_to_slots([records[part]
                                                            for part in parts])
                    for name, parts in assembly_plan.assemblies.items()
                }

        elif data.input_type == 'constructs':
            records = records_from_data_files(data.constructs_files)
            slots_dict = OrderedDict([
                (rec.id, kappagate.construct_record_to_slots(
                    rec, backbone_annotations=(data.backbone_string,)))
                for rec in records
            ])
       
        elif data.input_type == 'overhangs':
            overhangs = csv_to_list(data.overhangs_list)
            slots_dict = {
                'slots': kappagate.overhangs_list_to_slots(overhangs)}
        
        pdf_io = BytesIO()

        log(message="Running simulations...")
        summaries = []
        plain_summaries = {}
        for name, slots in log.iter_bar(construct=slots_dict.items()):
            proba, _, _ = kappagate.predict_assembly_accuracy(
                slots=slots,
                corrective_factor=data.dataset_corrective_factor,
                annealing_data=(data.dataset_tm, data.dataset_duration))
            plain_summaries[name] = kappagate.success_rate_facts(
                proba, plain_text=True)
            stats = kappagate.success_rate_facts(proba, plain_text=False)
            
            stats['construct_name'] = name
            summaries.append(OrderedDict([
                ('construct_name', name),
                ('good_clones_in_percent', stats['success_rate_percent']),
                ('pick_X_clones_for_95p_confidence', stats['min_trials_q95']),
                ('little_hope_after_X_bad_clones', stats['max_trials_q99'])
            ]))
        summary_table = pandas.DataFrame(summaries)
        summary_table.columns = list(summaries[0].values())
        csv_data = summary_table.to_csv()

        with PdfPages(pdf_io) as pdf:
            for name, slots in log.iter_bar(construct=slots_dict.items()):
                figsize = max(9, len(slots) - 4)
                fig, ax = plt.subplots(1, figsize=(figsize, figsize))
                kappagate.plot_circular_interactions(
                    slots, annealing_data=(data.dataset_tm,
                                           data.dataset_duration),
                    rate_limit=data.rate_limit_in_plots,
                    corrective_factor=data.dataset_corrective_factor,
                    ax=ax)
                ax.set_title(name, fontdict=dict(size=16, weight='bold'))
                wrapped = "\n".join(textwrap.wrap(plain_summaries[name], 80))
                ax.text(0, -0.1, wrapped, fontdict=dict(size=10),
                        transform=ax.transAxes)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        pdf_data = pdf_io.getvalue()
        html_pdf_data = data_to_html_data(pdf_data, 'pdf')

        return dict(
            csv_file_data=dict(
                name="predicted_valid_clone_rates.csv",
                data=csv_data, #data_to_html_data(csv_data.encode(), 'csv'),
                mimetype='application/csv'),
            pdf_data=html_pdf_data
        )

class PredictBadCloneRatesView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
