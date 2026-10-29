import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from complex_coherence.conservative import repair
from taxon_consistency.adjust_ortools import solve_sat_ortools_taxon_assignment

class ConservativeTests(unittest.TestCase):
    def test_no_evidence_abstains_and_preserves(self):
        scores = {'a': {'c': .9, 'root': .9}, 'b': {}}
        final, audit = repair(scores, {'c': {'root'}}, {'c'}, {})
        self.assertEqual(final, scores)
        self.assertEqual(audit['unresolved'][0]['reason'], 'no_authorized_edit')

    def test_evidenced_promotion_closes_ancestors(self):
        scores = {'a': {'c': .9, 'root': .9}, 'b': {}, 'z': {}}
        final, audit = repair(scores, {'c': {'root'}}, {'c'}, {('c','b'): 'independent evidence'})
        self.assertGreater(final['b']['c'], .001)
        self.assertGreater(final['b']['root'], .001)
        self.assertEqual(final['z'], {})
        self.assertEqual(audit['unresolved'], [])

    def test_protection_overrides_negative_conflict(self):
        with self.assertRaises(ValueError):
            repair({'a': {'c': .9}}, {}, {'c'}, {},
                   {('c', 'a'): 'negative'}, {('c', 'a')})

    def test_authorized_demotion_only(self):
        final, audit = repair({'a': {'c': .01}, 'b': {}}, {}, {'c'}, {},
                              {('c','a'): 'negative evidence'})
        self.assertEqual(audit['demotions'], [['c', 'a']])
        self.assertEqual(audit['promotions'], [])
        self.assertLessEqual(final['a']['c'], .001)

    def test_hierarchy_incompatible_candidate_abstains(self):
        final, audit = repair({'a': {'c': .9, 'root': .9}, 'b': {}},
                              {'c': {'root'}}, {'c'}, {('c','b'): 'support'},
                              singleton_penalty=.0001)
        self.assertEqual(final['b'], {})
        self.assertEqual(audit['unresolved'][0]['reason'], 'objective_or_coupled_constraints')

class UncertaintyTests(unittest.TestCase):
    def solve(self, costs, delta=0):
        result = {}
        solve_sat_ortools_taxon_assignment(
            {'p': {'g', 'h'}}, {'g':['NCBITaxon_2'], 'h':['NCBITaxon_3']}, {}, {}, costs,
            {'NCBITaxon_2': {'NCBITaxon_1'}, 'NCBITaxon_3': {'NCBITaxon_1'}}, {}, {},
            uncertainty=result, query_taxa={'NCBITaxon_2','NCBITaxon_3'}, near_optimal_units=delta)
        return result

    def test_alternative_go_repairs_are_ambiguous(self):
        r = self.solve({'g':1, 'h':1})
        self.assertEqual(r['required'], [])
        self.assertEqual(set(r['possible']), {'NCBITaxon_2','NCBITaxon_3'})
        self.assertEqual(r['optimum_units'], 1000000)

    def test_unique_optimum_and_near_optimal_alternative(self):
        r = self.solve({'g':2, 'h':1})
        self.assertEqual(r['required'], ['NCBITaxon_2'])
        self.assertEqual(self.solve({'g':2, 'h':1}, 1000000)['required'], [])


class RecoveryTests(unittest.TestCase):
    def test_recovery_inputs_do_not_use_hidden_membership(self):
        from analysis.revision_extensions import recovery_case
        original = {'seed': {'c':.9}, 'hidden': {'c':.8}, 'other': {'c':.1}}
        outputs, evidence = recovery_case(original, 'c', 'seed', .3, {}, {'c'}, 10)
        self.assertEqual(evidence, ['hidden'])
        self.assertEqual(outputs['conservative']['members'], ['hidden','seed'])
        self.assertTrue(outputs['demotion_only']['seed_removed'])

    def test_no_predictor_support_abstains(self):
        from analysis.revision_extensions import recovery_case
        outputs, evidence = recovery_case({'seed': {'c':.9}, 'hidden': {}}, 'c', 'seed', .3, {}, {'c'}, 10)
        self.assertEqual(evidence, [])
        self.assertTrue(outputs['conservative']['unresolved'])
        self.assertFalse(outputs['conservative']['seed_removed'])

class EdgeTests(unittest.TestCase):
    def test_reused_report_does_not_retain_certificates_after_unproven_solve(self):
        from unittest.mock import patch
        from ortools.sat.python import cp_model
        original = cp_model.CpSolver.Solve
        report = {'optimum_units': 0, 'checks': {'NCBITaxon_2': 'INFEASIBLE'},
                  'required': ['NCBITaxon_2']}
        def feasible(solver, model):
            original(solver, model)
            return cp_model.FEASIBLE
        with patch.object(cp_model.CpSolver, 'Solve', feasible):
            solve_sat_ortools_taxon_assignment(
                {'p': {'g'}}, {'g': ['NCBITaxon_2']}, {}, {}, {'g': 1},
                {'NCBITaxon_2': {'NCBITaxon_1'}}, {}, {},
                uncertainty=report, query_taxa={'NCBITaxon_2'})
        self.assertEqual(report['status'], 'FEASIBLE')
        self.assertEqual(report['required'], [])
        self.assertEqual(report['undetermined'], ['NCBITaxon_2'])
        self.assertNotIn('checks', report)
        self.assertNotIn('optimum_units', report)

    def test_low_level_protection_outside_domain_is_validated(self):
        from complex_coherence.adjust_ortools import solve_sat_ortools_optimized
        with self.assertRaises(ValueError):
            solve_sat_ortools_optimized({'a':set()}, set(), {}, {}, {},
                                       protected_pairs={('c','a')})

    def test_penalty_cannot_silently_round_to_zero(self):
        with self.assertRaises(ValueError):
            repair({'a': {'c':.9}}, {}, {'c'}, {}, singleton_penalty=1e-9)

    def test_no_complex_edits_preserves_unrelated_scores(self):
        scores = {'a': {'g':.2}, 'b':{}}
        result, audit = repair(scores, {}, {'c'}, {})
        self.assertEqual(result, scores)
        self.assertEqual(audit['flips'], 0)

    def test_unclosed_input_is_rejected(self):
        with self.assertRaises(ValueError):
            repair({'a': {'c':.9}}, {'c':{'root'}}, {'c'}, {})

    def test_query_timeout_is_not_entailment(self):
        from unittest.mock import patch
        from ortools.sat.python import cp_model
        original = cp_model.CpSolver.Solve
        calls = []
        def solve(solver, model):
            calls.append(1)
            return original(solver, model) if len(calls)==1 else cp_model.UNKNOWN
        with patch.object(cp_model.CpSolver, 'Solve', solve):
            result = UncertaintyTests().solve({'g':2, 'h':1})
        self.assertEqual(result['required'], [])
        self.assertEqual(set(result['undetermined']), {'NCBITaxon_2','NCBITaxon_3'})



class ExhaustiveTests(unittest.TestCase):
    def test_conservative_objective_matches_exhaustive_search(self):
        import itertools
        for candidate_score in [0, .1, .29]:
            for penalty in [.05, .4, 1.0]:
                scores = {'a': {'c':.9, 'root':.9},
                          'b': {'c':candidate_score, 'root':candidate_score}, 'z': {}}
                final, audit = repair(scores, {'c':{'root'}}, {'c'},
                                      {('c','b'):'support'}, threshold=.3,
                                      singleton_penalty=penalty)
                possible_costs = []
                # a is preserved, z has no evidence; b may add c and/or root.
                for child, parent in itertools.product([0,1], repeat=2):
                    if child > parent:
                        continue
                    flip_cost = max(1, int(abs(candidate_score-.3)*1000000))
                    possible_costs.append((child+parent)*flip_cost +
                                          (1-child)*round(penalty*1000000))
                self.assertEqual(round(audit['solver']['objective']*1000000), min(possible_costs))



class CorroborationTests(unittest.TestCase):
    def test_corroboration_can_recover_an_unscored_candidate(self):
        from analysis.revision_extensions import recovery_case
        outputs, _ = recovery_case({'seed': {'c':.9}, 'hidden': {}, 'other': {}},
                                    'c', 'seed', .3, {}, {'c'}, 10,
                                    {('c','hidden'): 'both other predictors'})
        self.assertEqual(outputs['conservative']['members'], ['seed'])
        self.assertEqual(outputs['conservative_corroborated']['members'], ['hidden','seed'])

    def test_other_predictors_must_both_support_same_pair(self):
        import tempfile
        from types import SimpleNamespace
        from analysis.revision_extensions import corroborating_evidence
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for predictor, content in [('seq_sim_timeset','p\tc|0.01\nq\tc|0.01\n'),
                                       ('sprof-go','p\tc|0.001\nq\tc|0.02\n')]:
                path = root/predictor/'cc'/'stage1'
                path.mkdir(parents=True)
                (path/'optimized_fold_01_taxon_1.tsv').write_text(content)
            evidence, sources = corroborating_evidence(SimpleNamespace(
                predictor='deepgo-se', taxon='1', results=root), {'c'})
            self.assertEqual(set(evidence), {('c','q')})
            self.assertEqual(len(sources), 2)

if __name__ == '__main__':
    unittest.main()
