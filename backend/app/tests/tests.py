"""Tests for the backend"""

from .tools import logprint, AppTestCase, load_file_to_dict

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


class SimulateGGAssemblies(AppTestCase):
    endpoint = 'simulate_gg_assemblies'
    defaults = dict(
        enzyme='Autoselect',
        parts=[],
        connectors=[],
        select_connectors=False,
        include_fragments=False,
        use_assembly_plan=False,
        single_assemblies=False,
        assembly_plan=None
    )

    def test_dual_assembly(self):
        parts = [
            load_file_to_dict(['simulate_gg_assemblies', "%s.gb" % p])
            for p in ['partA', 'partA2', 'partB', 'partC', 'receptor']
        ]
        for part in parts:
            part['circularity'] = True
        logprint(parts)
        response = self.run_job(parts=parts)#parts=parts)
        logprint(response)
        self.assertEqual(response['infos']['nconstructs'], 2)
