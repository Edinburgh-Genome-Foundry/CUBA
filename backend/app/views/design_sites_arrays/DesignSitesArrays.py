"""Bla."""

from rest_framework import serializers
import flametree
from zymp import (stacked_sites_array, plot_sequence_sites,
                  annotate_enzymes_sites, write_record)
from ..base import AsyncWorker, StartJobView
from ..tools import data_to_html_data, file_to_filelike_object, csv_to_list
from ..serializers import FileSerializer


class serializer_class(serializers.Serializer):
    """Serializer."""
    sites_to_include = serializers.CharField()
    forbidden_sites = serializers.CharField(allow_blank=True)
    mandatory_sites = serializers.CharField(allow_blank=True)
    enforce_unique_sites = serializers.BooleanField()

class worker_class(AsyncWorker):

    def work(self):
        data = self.data
        self.logger(message='analyzing the data...')

        sites_to_include = csv_to_list(data.sites_to_include)
        forbidden_sites= csv_to_list(data.forbidden_sites)
        mandatory_sites = set(csv_to_list(data.mandatory_sites))

        def success_condition(seq, inseq, notinseq):
            return mandatory_sites <= set(inseq)

        # DESIGN AN OPTIMIZED SEQUENCE WITH ZYMP
        seq, sites_in_seq, leftover = stacked_sites_array(
                sites_to_include, forbidden_enzymes=forbidden_sites,
                unique_sites=data.enforce_unique_sites,
                success_condition=success_condition, tries=200)
        
        zip_root = flametree.file_tree("@memory")

        # PLOT A SUMMARY
        ax = plot_sequence_sites(seq, sites_in_seq)
        ax.figure.savefig(zip_root._file("stacked_array.pdf").open('wb'),
                        format='pdf', bbox_inches='tight')

        # WRITE THE SEQUENCE AND SITE ANNOTATIONS AS A RECORD
        record = annotate_enzymes_sites(
            seq, sites_in_seq, forbidden_enzymes=forbidden_sites)
        write_record(record, zip_root._file("record.gb").open('w'))
        message = ("Generated a %dn sequence with %d sites."
                   % (len(seq), len(sites_in_seq)))
        if len(leftover):
            message += "Sites non included: %s" % ", ".join(leftover)
        figure_data = zip_root["stacked_array.pdf"].read("rb")
        return dict(
            success=True,
            message=message,
            genbank_file= dict(
                data=data_to_html_data(zip_root["record.gb"].read('r').encode(),
                                       'genbank'),
                name='stacked_site_array.gb',
                mimetype='application/genbank'
            ),
            figure_data=data_to_html_data(figure_data, 'pdf')
    )
        

class DesignSitesArraysView(StartJobView):
    serializer_class = serializer_class
    worker_class = worker_class
