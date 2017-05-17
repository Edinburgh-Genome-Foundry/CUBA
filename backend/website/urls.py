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
    url(r'^poll$', views.PollJobView.as_view()),

    url(r'^start/predict_digests$',
        views.PredictDigestsView.as_view()),
    url(r'^start/simulate_cloning$',
        views.SimulateCloningView.as_view()),
    url(r'^start/select_digestions$',
        views.SelectDigestionsView.as_view()),
    url(r'^start/evaluate_manufacturability$',
        views.EvaluateManufacturabilityView.as_view()),

    url(r'^docs/', include('rest_framework_docs.urls')),
    url(r'^django-rq/', include('django_rq.urls')),
    url(r'^admin/', admin.site.urls)
]
