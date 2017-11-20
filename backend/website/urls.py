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
    url(r'^api/poll$', views.PollJobView.as_view()),

    url(r'^api/start/analyze_digests$',
        views.AnalyzeDigestsView.as_view()),
    url(r'^api/start/convert_sequence_files$',
        views.ConvertSequenceFilesView.as_view()),
    url(r'^api/start/design_overhangs$',
        views.DesignOverhangsView.as_view()),
    url(r'^api/start/evaluate_manufacturability$',
        views.EvaluateManufacturabilityView.as_view()),
    url(r'^api/start/find_common_blocks$',
        views.FindCommonBlocksView.as_view()),
    url(r'^api/start/predict_digestions$',
        views.PredictDigestionsView.as_view()),
    url(r'^api/start/sculpt_a_sequence$',
        views.SculptASequenceView.as_view()),
    url(r'^api/start/select_digestions$',
        views.SelectDigestionsView.as_view()),
    url(r'^api/start/select_primers$',
        views.SelectPrimersView.as_view()),
    url(r'^api/start/simulate_gg_assemblies$',
        views.SimulateGGAssembliesView.as_view()),
    url(r'^api/start/sketch_constructs$',
        views.SketchConstructsView.as_view()),
    url(r'^api/start/swap_donor_vector_part$',
        views.SwapDonorVectorPartView.as_view()),

    url(r'^api/docs/', include('rest_framework_docs.urls')),
    url(r'^api/django-rq/', include('django_rq.urls')),
    url(r'^api/admin/', admin.site.urls)
]
