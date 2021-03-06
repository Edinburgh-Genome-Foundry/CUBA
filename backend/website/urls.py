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
    url(r'^api/start/compare_two_sequences$',
        views.CompareTwoSequencesView.as_view(),
        name='analyze_digests'),
    url(r'^api/start/convert_sequence_files$',
        views.ConvertSequenceFilesView.as_view(),
        name='convert_sequence_files'),
    url(r'^api/start/create_assembly_picklists$',
        views.CreateAssemblyPicklistsView.as_view(),
        name='create_assembly_picklists'),
    url(r'^api/start/design_part_test_batches$',
        views.DesignPartTestBatchesView.as_view(),
        name='design_part_test_batches'),
    url(r'^api/start/design_sites_arrays$',
        views.DesignSitesArraysView.as_view(),
        name='design_sites_arrays'),
    url(r'^api/start/design_overhangs$',
        views.DesignOverhangsView.as_view(),
        name='design_overhangs'),
    url(r'^api/start/domesticate_part_batches$',
        views.DomesticatePartBatchesView.as_view(),
        name='domesticate_part_batches'),
    url(r'^api/start/evaluate_manufacturability$',
        views.EvaluateManufacturabilityView.as_view(),
        name='evaluate_manufacturability'),
    url(r'^api/start/find_common_blocks$',
        views.FindCommonBlocksView.as_view(),
        name='find_common_blocks'),
    url(r'^api/start/insert_parts_on_backbones$',
        views.InsertPartsOnBackbonesView.as_view(),
        name='insert_parts_on_backbones'),
    url(r'^api/start/find_saboteur_parts$',
        views.FindSaboteurPartsView.as_view(),
        name='find_saboteur_parts'),
    url(r'^api/start/plot_sequence_features$',
        views.PlotSequenceFeaturesView.as_view(),
        name='plot_sequence_features'),
    url(r'^api/start/predict_digestions$',
        views.PredictDigestionsView.as_view(),
        name='predict_digestions'),
    url(r'^api/start/predict_bad_clone_rates$',
        views.PredictBadCloneRatesView.as_view(),
        name='predict_digestions'),
    url(r'^api/start/rearray_plates$',
        views.RearrayPlatesView.as_view(),
        name='rearray_plates'),
    url(r'^api/start/render_sequenticons$',
        views.RenderSequenticonsView.as_view(),
        name='render_sequenticons'),
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
    url(r'^api/start/simulate_multi_method_assemblies$',
        views.SimulateMultiMethodAssembliesView.as_view(),
        name='simulate_multi_method_assemblies'),
    url(r'^api/start/sketch_constructs$',
        views.SketchConstructsView.as_view(),
        name='sketch_constructs'),
    url(r'^api/start/substitute_parts_overhangs$',
        views.SubstitutePartsOverhangsView.as_view(),
        name='sketch_constructs'),
    url(r'^api/start/transfer_features$',
        views.TransferFeaturesView.as_view(),
        name='transfer_features'),
    url(r'^api/start/view_overhangs_crosstalk$',
        views.ViewOverhangsCrosstalkView.as_view(),
        name='view_overhangs_crosstalk'),
    url(r'^api/start/extract_features$',
        views.ExtractFeaturesView.as_view(),
        name='extract_features'),

    url(r'^api/docs/', include('rest_framework_docs.urls')),
    url(r'^api/django-rq/', include('django_rq.urls')),
    url(r'^api/admin/', admin.site.urls)
]
