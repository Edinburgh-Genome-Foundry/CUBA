"""Bla."""

from rest_framework import serializers
import kappagate
from ..serializers import FileSerializer
from ..base import AsyncWorker, StartJobView
from ..tools import (
    matplotlib_figure_to_svg_base64_data,
    records_from_data_files,
)
import tatapov


class serializer_class(serializers.Serializer):
    """Serializer."""

    overhangs = serializers.CharField(allow_blank=True)
    construct = FileSerializer(allow_null=True)
    temperature = serializers.CharField()
    incubation = serializers.CharField()


class worker_class(AsyncWorker):
    def work(self):
        data = self.data
        overhangs = [
            s.strip().upper()
            for s in data.overhangs.split(",")
            if len(s.strip())
        ]
        if overhangs == []:
            record = records_from_data_files([data.construct])[0].upper()
            slots = kappagate.construct_record_to_slots(record, backbone_annotations=("",))
            overhangs = [o for s in slots for o in s if set(o) <= set("ATGC")]

        temperature, incubation = data.temperature, data.incubation
        data = tatapov.annealing_data[temperature][incubation]
        subset = tatapov.data_subset(data, overhangs, add_reverse=True)
        ax, _ = tatapov.plot_data(subset, figwidth=2 + 0.5 * len(overhangs))
        ax.figure.tight_layout()
        fig_data = matplotlib_figure_to_svg_base64_data(
            ax.figure, bbox_inches="tight"
        )
        return {"figure_data": fig_data, "success": True}


class ViewOverhangsCrosstalkView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
