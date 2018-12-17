"""Bla."""

from rest_framework import serializers
from ..base import AsyncWorker, StartJobView
from ..tools import matplotlib_figure_to_svg_base64_data
import tatapov

class serializer_class(serializers.Serializer):
    """Serializer."""
    overhangs = serializers.CharField()
    temperature = serializers.CharField()
    incubation = serializers.CharField()

class worker_class(AsyncWorker):

    def work(self):
        overhangs = [
            s.strip().upper()
            for s in self.data.overhangs.split(',')
            if len(s.strip())
        ]
        temperature, incubation = self.data.temperature, self.data.incubation
        data = tatapov.annealing_data[temperature][incubation]
        subset = tatapov.data_subset(data, overhangs, add_reverse=True)
        ax, _ = tatapov.plot_data(subset, figwidth=2 + 0.5 * len(overhangs))
        ax.figure.tight_layout()
        fig_data = matplotlib_figure_to_svg_base64_data(ax.figure,
                                                        bbox_inches="tight")
        return {
          'figure_data': fig_data,
          'success': True
        }

class ViewOverhangsCrosstalkView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
