"""Regression cases for the reviewer-revision audit; no external data needed."""
import sys
import unittest
import tempfile
import os
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from complex_coherence import adjust_ortools as cc
from taxon_consistency import adjust_ortools as tc


class ScoreContractTests(unittest.TestCase):
    def test_sparse_linear_flip_objective_matches_exhaustive_optimum(self):
        import itertools
        ann={'p1':{'complex','root'},'p2':{'root'}}
        hierarchy={'complex':{'root'}}
        pairs=sorted((g,p) for g in ['complex','root'] for p in ann)
        costs={'p1':{'complex':.6,'root':.7},'p2':{'complex':.01,'root':.8}}
        for demotion_only, penalty in [(False,None),(True,None),(False,.005),(False,.2)]:
            candidates=[]
            for bits in itertools.product([False,True],repeat=len(pairs)):
                values=dict(zip(pairs,bits))
                if any(values['complex',p] and not values['root',p] for p in ann): continue
                if demotion_only and values['complex','p2']:continue
                singleton=sum(values['complex',p] for p in ann)==1
                if penalty is None and singleton:continue
                objective=sum(max(1,int(costs[p][g]*1e6)) for (g,p),v in values.items() if v!=(g in ann[p]))
                if penalty is not None:objective+=round(penalty*1e6)*singleton
                candidates.append(objective/1e6)
            stats={}
            result,_=cc.solve_sat_ortools_optimized(ann,set(pairs),hierarchy,costs,{'complex':['p1','p2']},
                demotion_only=demotion_only,violation_penalty=penalty,solver_stats=stats)
            self.assertEqual(stats['status'],'OPTIMAL')
            self.assertAlmostEqual(stats['objective'],min(candidates))

    def test_curated_location_is_not_counted_as_complex_membership(self):
        import json
        import io
        from contextlib import redirect_stdout
        from unittest.mock import patch
        from analysis.complex_validation import validate_stage2 as validator
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as directory:
            base=Path(directory); before=base/'before'; after=base/'after'; before.mkdir(); after.mkdir()
            name='optimized_fold_0_taxon_1.tsv'
            (before/name).write_text('P12345\n')
            (after/name).write_text('P12345\tGO:0000001|0.301\tGO:0005737|0.301\n')
            index=base/'index.json'; index.write_text(json.dumps({'1::GO:0000001':['P12345'],'1::GO:0005737':['P12345']}))
            go=base/'go.obo'
            go.write_text('format-version: 1.2\n\n[Term]\nid: GO:0032991\nname: complex\nnamespace: cellular_component\n\n[Term]\nid: GO:0000001\nname: child complex\nnamespace: cellular_component\nis_a: GO:0032991 ! complex\n\n[Term]\nid: GO:0005737\nname: cytoplasm\nnamespace: cellular_component\n')
            exclusions=base/'exclude.tsv'; exclusions.write_text('GO_term\tClassification\n')
            output=base/'result.json'
            argv=['validate','--stage1_dir',str(before),'--stage2_dir',str(after),'--subunit_index',str(index),
                '--go',str(go),'--complex-exclusions',str(exclusions),'--threshold','.3','--output_file',str(output)]
            with patch.object(sys,'argv',argv), redirect_stdout(io.StringIO()):
                validator.main()
            report=json.loads(output.read_text())
            self.assertEqual(report['promotions']['total'],2)
            self.assertEqual(report['promotions']['with_curated_reference'],1)
            self.assertEqual(report['promotions']['true_subunit_promotions'],1)

    def test_public_pipeline_reports_empty_active_denominator(self):
        from unittest.mock import patch
        from pipeline import run_adjustment_pipeline as pipeline
        from contextlib import redirect_stdout
        import io
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as directory:
            base=Path(directory)
            (base/'predictions_fold_0_taxon_1.tsv').write_text('p\n')
            args=dict(predictions_dir=str(base),optimized_dir=str(base/'out'),fold_thresholds={'0':.3},
                constraints_file='unused',go_hierarchy_file='unused',taxon_hierarchy_file='unused',
                ncbitaxon_hierarchy_file='unused',complexes_file='unused')
            with patch.object(pipeline,'load_constraints',return_value=({},{})), \
                 patch.object(pipeline,'load_go_hierarchy',return_value={}), \
                 patch.object(pipeline,'load_taxon_hierarchy',return_value=({},{},{})), \
                 patch.object(pipeline,'adjust_per_taxon',return_value=(0,0)):
                output=io.StringIO()
                with redirect_stdout(output):
                    self.assertTrue(pipeline.run_adjustment_pipeline(**args))
                self.assertIn('undefined (no active inputs)',output.getvalue())
                with self.assertRaises(ValueError):
                    pipeline.run_adjustment_pipeline(**args,skip_existing=True)

    def test_export_validator_checks_scores_and_metadata_counts(self):
        from analysis.validate_repair_exports import check_pair
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as directory:
            before, after = Path(directory)/'before.tsv', Path(directory)/'after.tsv'
            before.write_text('p1\tcomplex|0.9\troot|0.9\np2\troot|0.9\n')
            after.write_text('p1\tcomplex|0.9\troot|0.9\np2\tcomplex|0.301\troot|0.9\n')
            result=check_pair(before,after,.3,{'complex': {'root'}},{'complex'},set(),False)
            self.assertEqual(result['promotions'],1)
            self.assertEqual(result['demotions'],0)
            with self.assertRaises(ValueError):
                check_pair(before,after,.3,{'complex': {'root'}},{'complex'},set(),True)
            after.write_text('p1\tcomplex|0.9\troot|0.9\np2\tcomplex|0.3\troot|0.9\n')
            with self.assertRaises(ValueError):
                check_pair(before,after,.3,{'complex': {'root'}},{'complex'},set(),False)

    def test_soft_taxon_penalty_is_optimized(self):
        isa = defaultdict(set, {'NCBITaxon_2': {'NCBITaxon_1'}, 'NCBITaxon_3': {'NCBITaxon_1'}})
        for penalty, expected in [(0, True), (.2, True), (.9, False)]:
            result, _ = tc.solve_sat_ortools_hierarchy(
                {'p': {'g'}}, {'g': ['NCBITaxon_2']}, {}, defaultdict(set),
                {'g': .7}, isa, {}, {}, '3', violation_costs={'g': penalty})
            self.assertEqual(result.get('p', {}).get('g', True), expected)

    def test_union_members_are_not_assumed_disjoint(self):
        from ortools.sat.python import cp_model
        model = cp_model.CpModel()
        isa = defaultdict(set, {'NCBITaxon_33208': {'NCBITaxon_2759'}, 'NCBITaxon_2759': {'NCBITaxon_1'}})
        union = {'NCBITaxon_Union_1': {'NCBITaxon_33208', 'NCBITaxon_2759'}}
        variables = tc._add_taxon_model(model, {'NCBITaxon_Union_1'}, isa, {}, union, '33208')
        solver = cp_model.CpSolver()
        self.assertEqual(solver.Solve(model), cp_model.OPTIMAL)
        self.assertTrue(solver.BooleanValue(variables['NCBITaxon_Union_1']))
        self.assertTrue(solver.BooleanValue(variables['NCBITaxon_33208']))
        self.assertTrue(solver.BooleanValue(variables['NCBITaxon_2759']))

    def test_soft_complex_penalty_is_optimized_not_applied_greedily(self):
        ann = {'p': {'complex'}}
        pairs = {('complex', 'p')}
        costs = {'p': {'complex': .7}}
        for penalty, expected in [(.2, True), (.9, False)]:
            result, _ = cc.solve_sat_ortools_optimized(
                ann, pairs, defaultdict(set), costs, {'complex': ['p']},
                violation_penalty=penalty)
            self.assertEqual(result['p']['complex'], expected)

    def test_uniprot_entry_name_requires_a_mapping(self):
        from analysis.complex_validation.validate_stage2 import protein_to_accession
        self.assertIsNone(protein_to_accession('143BA_DANRE'))
        self.assertEqual(protein_to_accession('143BA_DANRE', {'143BA_DANRE': 'Q9PUL6'}), 'Q9PUL6')
        self.assertEqual(protein_to_accession('A0A6P5QAZ4_MUSCR'), 'A0A6P5QAZ4')

    def test_stage1_removes_genome_term_and_descendants_at_all_scores(self):
        raw = {'p1': {'bad': .9, 'ok': .8}, 'p2': {'bad': .1, 'child': .01}}
        out = tc.get_final_predictions({'p1': {'bad': False}}, raw, {'child': {'bad'}})
        self.assertEqual(out, {'p1': {'ok': .8}, 'p2': {}})

    def test_stage2_preserves_scores_outside_active_vocabulary(self):
        raw = {'p1': {'complex': .9, 'low': .02}, 'p2': {'complex': .29}}
        out = cc.get_final_predictions({'p1': {'complex': False}}, raw, {'complex'}, .3)
        self.assertEqual(out['p1'], {'complex': .3, 'low': .02})
        self.assertEqual(out['p2'], raw['p2'])

    def test_active_score_just_above_threshold_is_unchanged(self):
        raw = {'p': {'g': .30004}}
        self.assertEqual(cc.get_final_predictions({'p': {'g': True}}, raw, {'g'}, .3), raw)

    def test_isolated_complex_is_constrained_in_full_model(self):
        result, _ = cc.solve_sat_ortools_hierarchy({'p': {'complex'}}, {}, {'p': {'complex': .8}}, {'complex'})
        self.assertFalse(result['p']['complex'])

    def test_candidate_threshold_keeps_the_singleton_available(self):
        self.assertEqual(cc.get_participating_proteins({'p': {'complex'}}, {'complex'},
                         {'p': {'complex': .9}}, participating_threshold=.95), {'complex': ['p']})

    def test_promoted_ancestor_is_exported(self):
        out = cc.get_final_predictions({'p1': {'parent': True}}, {'p1': {'child': .9}}, {'child'}, .3)
        self.assertGreater(out['p1']['parent'], .3)

    def test_demotion_only_is_enforced_by_both_solvers(self):
        hierarchy = defaultdict(set, {'complex': {'root'}})
        raw = {'p1': {'root': .9, 'complex': .9}, 'p2': {'root': .9, 'complex': .29}}
        active, terms = cc.get_annotated_predictions(raw, .3)
        costs, _ = cc.compute_flip_cost(raw, .3)
        full, _ = cc.solve_sat_ortools_hierarchy(active, hierarchy, costs, {'complex'}, demotion_only=True)
        participants = cc.get_participating_proteins(active, {'complex'}, raw, demotion_only=True)
        pairs = cc.get_go_protein_pairs(active, hierarchy, participants)
        sparse, _ = cc.solve_sat_ortools_optimized(active, pairs, hierarchy, costs, participants, demotion_only=True)
        for adjusted in [full, sparse]:
            out = cc.get_final_predictions(adjusted, raw, terms, .3)
            self.assertLessEqual(out['p1']['complex'], .3)
            self.assertEqual(out['p2']['complex'], .29)

    def test_one_protein_cannot_support_a_heteromeric_complex(self):
        active = {'p1': {'complex', 'root'}}
        costs = defaultdict(lambda: defaultdict(lambda: .3))
        out, _ = cc.solve_sat_ortools_hierarchy(active, defaultdict(set, {'complex': {'root'}}), costs, {'complex'})
        self.assertFalse(out['p1']['complex'])


if __name__ == '__main__':
    unittest.main()
