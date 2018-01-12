"""website URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""

import matplotlib
matplotlib.use('Agg')

from django.conf.urls import url, include
from django.contrib import admin
import app.views as views

urlpatterns = [
    url(r'^api/poll$', views.PollJobView.as_view(), name='poll'),

    url(r'^api/start/analyze_digests$',
        views.AnalyzeDigestsView.as_view(),
        name='analyze_digests'),
    url(r'^api/start/convert_sequence_files$',
        views.ConvertSequenceFilesView.as_view(),
        name='convert_sequence_files'),
    url(r'^api/start/design_overhangs$',
        views.DesignOverhangsView.as_view(),
        name='design_overhangs'),
    url(r'^api/start/evaluate_manufacturability$',
        views.EvaluateManufacturabilityView.as_view(),
        name='evaluate_manufacturability'),
    url(r'^api/start/find_common_blocks$',
        views.FindCommonBlocksView.as_view(),
        name='find_common_blocks'),
    url(r'^api/start/predict_digestions$',
        views.PredictDigestionsView.as_view(),
        name='predict_digestions'),
    url(r'^api/start/sculpt_a_sequence$',
        views.SculptASequenceView.as_view(),
        name='sculpt_a_sequence'),
    url(r'^api/start/select_digestions$',
        views.SelectDigestionsView.as_view(),
        name='select_digestions'),
    url(r'^api/start/select_primers$',
        views.SelectPrimersView.as_view(),
        name='select_primers'),
    url(r'^api/start/simulate_gg_assemblies$',
        views.SimulateGGAssembliesView.as_view(),
        name='simulate_gg_assemblies'),
    url(r'^api/start/sketch_constructs$',
        views.SketchConstructsView.as_view(),
        name='sketch_constructs'),
    url(r'^api/start/swap_donor_vector_part$',
        views.SwapDonorVectorPartView.as_view(),
        name='swap_donor_vector_part'),

    url(r'^api/docs/', include('rest_framework_docs.urls')),
    url(r'^api/django-rq/', include('django_rq.urls')),
    url(r'^api/admin/', admin.site.urls)
]
