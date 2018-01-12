
"""
'sequence': 'None',
'goal': 'overhangs_set',
'overhangs_differences': '2',
'gc_content': '[25, 75]',
'mandatory_overhangs': '',
'forbidden_overhangs': '',
'cutting_mode': 'equal',
'extremities': 'True',
'n_overhangs': '10',
'n_fragments': '2',
'auto_overhangs': 'True',
'left_flank_sequence': '',
'right_flank_sequence': '',
'allow_edits': 'False',
'domain_name': 'django:8082'
"""

from .tools import logprint, AppTestCase

class DesignOverhangsTests(AppTestCase):
    endpoint = 'design_overhangs'
    defaults = dict(
        gc_content=[25, 75],
        overhangs_differences=2,
        mandatory_overhangs='',
        forbidden_overhangs='',
        cutting_mode='equal',
        auto_overhangs=False,
        extremities=True,
        n_overhangs=10,
        n_fragments=0,
        allow_edits=False,
        left_flank_sequence='',
        right_flank_sequence=''
    )

    def test_create_10_collection(self):
        response = self.run_job(goal='overhangs_set', n_overhangs=10,
                                auto_overhangs=False)
        self.assertEqual(len(response['overhangs']), 10)
