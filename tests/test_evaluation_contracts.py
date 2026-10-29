"""Scientific evaluation regression tests; uses a tiny synthetic ontology."""
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.evaluation_scope import select_annotations
from pipeline.genome_scores_evaluator import Ontology, evaluate_prediction_metrics

OBO = '''format-version: 1.2

[Term]
id: GO:0005575
name: cellular component
namespace: cellular_component

[Term]
id: GO:0000001
alt_id: GO:9000001
name: true label
namespace: cellular_component
is_a: GO:0005575

[Term]
id: GO:0000002
name: false label
namespace: cellular_component
is_a: GO:0005575
relationship: regulates GO:0000001

[Term]
id: GO:0000003
name: part
namespace: cellular_component
relationship: part_of GO:0000001

[Term]
id: GO:9000002
alt_id: GO:9000003
name: obsolete label
namespace: cellular_component
is_obsolete: true
'''


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home() / 'tmp')))
        path = Path(self.temp.name) / 'go.obo'
        path.write_text(OBO)
        self.go_file = path
        self.go = Ontology(path)

    def tearDown(self):
        self.temp.cleanup()

    def test_alternative_ids_merge_without_creating_false_positives(self):
        scores={'p': {'GO:9000001':1.0,'GO:0000001':.8}}
        result=evaluate_prediction_metrics(self.go,'cc',{'p':{'GO:9000001'}},scores)
        self.assertEqual(result['fmax'],1)
        self.assertEqual(result['aupr'],1)
        self.assertEqual(self.go.get_ancestors('GO:9000001'),{'GO:0000001','GO:0005575'})
        self.assertNotIn('GO:9000001',self.go.ont)
        audit={}
        scores['p']['GO:9000003']=.5
        self.assertEqual(self.go.normalize_scores(scores,audit=audit),{'p':{'GO:0000001':1.0}})
        self.assertEqual(audit['alternative_id_pairs'],1)
        self.assertEqual(audit['merged_duplicate_pairs'],1)
        self.assertEqual(audit['unknown_or_obsolete_pairs'],1)
        self.assertEqual(self.go.normalize_constraints({'GO:9000001':{'a'},'GO:0000001':{'b'}}),{'GO:0000001':{'a','b'}})

    def test_obsolete_aliases_and_invalid_scores_are_rejected(self):
        self.assertIsNone(self.go.canonical_id('GO:9000003'))
        self.assertEqual(self.go.get_ancestors('GO:9000003'),set())
        for score in [float('nan'),float('inf'),-.1,1.1]:
            with self.assertRaises(ValueError):
                self.go.normalize_scores({'p':{'GO:0000001':score}})

    def test_common_baseline_materializes_ancestors_without_demotion(self):
        raw=self.go.normalize_scores({'p':{'GO:9000001':.9,'GO:0005575':.2},'empty':{}})
        audit={}; baseline=self.go.propagate_scores(raw,audit)
        self.assertEqual(baseline,{'p':{'GO:0000001':.9,'GO:0005575':.9},'empty':{}})
        self.assertEqual(audit,{'hierarchy_added_pairs':0,'hierarchy_raised_pairs':1})
        baseline=self.go.propagate_scores({'p':{'GO:0000003':.8}},audit)
        self.assertEqual(baseline['p'],{'GO:0000003':.8,'GO:0000001':.8,'GO:0005575':.8})
        self.assertEqual(audit['hierarchy_added_pairs'],2)
        self.assertEqual(self.go.propagate_scores(baseline),baseline)

    def test_structural_report_is_saved_when_no_test_labels_match(self):
        from pipeline.genome_scores_evaluator import main
        output = Path(self.temp.name)/'report.json'
        structural = {'satisfiable': True, 'complete_has_part_percentage': 0, 'complex_coherence': 100}
        with patch('pipeline.genome_scores_evaluator.GAEF_evaluation', return_value=structural), patch('pipeline.genome_scores_evaluator.load_annotations', return_value={}):
            main('42', None, 'unused.tsv', 'unused', self.go_file, subontology='cc', output_file=str(output))
        report = json.loads(output.read_text())
        self.assertEqual(report['complex_coherence'], 100)
        self.assertEqual(report['prediction_metrics'], {})
        self.assertEqual(report['n_reference_proteins'], 0)

    def test_specific_terms_do_not_depend_on_iteration_order(self):
        from pipeline.genome_scores_evaluator import get_specific_with_ontology
        class OrderedSet(set):
            def __init__(self, values):
                super().__init__(values)
                self.order = values
            def __iter__(self):
                return iter(self.order)
        for order in [('GO:0000003', 'GO:0000001'), ('GO:0000001', 'GO:0000003')]:
            result = get_specific_with_ontology({'p': OrderedSet(order)}, self.go)
            self.assertEqual(result['p'], {'GO:0000003'})

    def test_score_delta_roundtrip_preserves_zero_and_missing_values(self):
        from analysis.aggregate_revision_audit import score_delta, apply_score_delta
        base = {'a': .8, 'b': .1, 'zero': 0.0}
        edited = {'a': .3, 'new': .301, 'zero': 0.0}
        self.assertEqual(apply_score_delta(base, score_delta(base, edited)), edited)
        self.assertEqual(apply_score_delta(base, {}), base)

    def test_relation_semantics(self):
        self.assertNotIn('GO:0000001', self.go.get_ancestors('GO:0000002'))
        self.assertIn('GO:0000001', self.go.get_ancestors('GO:0000003'))

    def test_reference_builder_does_not_inherit_regulated_function(self):
        from data_prep import obtain_annotations as builder
        ontology = builder.Ontology(self.go.filename, with_rels=True)
        labels = builder.propagate_annotations([['GO:0000002'], ['GO:0000003']], ontology)
        self.assertNotIn('GO:0000001', labels[0])
        self.assertIn('GO:0000001', labels[1])

    def test_reference_rebuild_records_source_and_removes_overpropagation(self):
        import argparse
        import gzip
        from analysis.rebuild_reference_annotations import rebuild
        root = Path(self.temp.name)
        source = root/'sources'/'Eukaryota'/'UP_TEST'
        source.mkdir(parents=True)
        with gzip.open(source/'UP_TEST_1.dat.gz', 'wt') as fh:
            fh.write('ID   P_TEST\nAC   Q12345;\nDR   GO; GO:0000002; C:test; EXP:UniProtKB.\n//\n')
        original = root/'original'; original.mkdir()
        (original/'annots_taxon_1.tsv').write_text('P_TEST\tGO:0000002\tGO:0000001\tGO:0005575\n')
        index = root/'index.tsv'; index.write_text('TaxonID\tProteomeID\n1\tUP_TEST\n')
        out = root/'rebuilt'
        manifest = rebuild(argparse.Namespace(source_root=root/'sources', original_annotations=original,
            proteome_index=index, go=self.go_file, output=out))
        row = manifest['records'][0]
        self.assertEqual(row['old_extra_pairs'], 1)
        self.assertEqual(row['new_reference_proteins_not_in_original_rows'], [])
        self.assertNotIn('GO:0000001', (out/'annots_taxon_1.tsv').read_text())
        self.assertEqual(len(row['source_sha256']), 64)
        self.assertTrue((out/'reference-manifest.json').exists())

    def test_reference_reuse_requires_score_coverage_in_every_condition(self):
        from analysis.aggregate_revision_audit import checked_reference_labels
        root = Path(self.temp.name)
        filename = 'annots_taxon_1.tsv'
        (root/filename).write_text('p\tGO:0000001\n')
        reference = dict(reference_file=filename, original_sha256='original', reference_sha256='rebuilt')
        data = dict(taxon='1', predictor='mlp_heldout', aspect='cc',
                    inputs={filename:'original'}, runs=[dict(test_scores={'p':{}})])
        self.assertEqual(set(checked_reference_labels(data, root, reference, {})), {'p'})
        data['runs'].append(dict(test_scores={}))
        with self.assertRaisesRegex(ValueError, 'Rerun required'):
            checked_reference_labels(data, root, reference, {})
        data['inputs'][filename] = 'unidentified'
        with self.assertRaisesRegex(ValueError, 'Unidentified'):
            checked_reference_labels(data, root, reference, {})

    def test_false_positive_not_in_any_test_label_counts(self):
        annotations = {'p': {'GO:0000001', 'GO:0005575'}}
        good = evaluate_prediction_metrics(self.go, 'cc', annotations, {'p': {'GO:0000001': .8}})
        bad = evaluate_prediction_metrics(self.go, 'cc', annotations, {'p': {'GO:0000001': .8, 'GO:0000002': .9}})
        self.assertEqual(good['fmax'], 1)
        self.assertAlmostEqual(bad['fmax'], 2/3)

    def test_perfect_unit_scores_have_unit_aupr(self):
        result = evaluate_prediction_metrics(self.go, 'cc', {'p': {'GO:0000001'}},
                                             {'p': {'GO:0000001': 1.0}})
        self.assertEqual(result['aupr'], 1)
        self.assertIsNone(result['avg_auc'])  # One protein cannot define ROC AUC.

    def test_temporal_mask_is_exact(self):
        annotations = {'new': {'a'}, 'old': {'b'}}
        self.assertEqual(select_annotations(annotations, '42', {'42': {'new'}}), {'new': {'a'}})
        self.assertEqual(select_annotations(annotations, '43', {'42': {'new'}}), {})

    def test_absent_scores_are_not_predictions_at_zero(self):
        result = evaluate_prediction_metrics(self.go, 'cc', {'p': {'GO:0000001', 'GO:0005575'}}, {})
        self.assertEqual(result['fmax'], 0)
        self.assertEqual(result['aupr'], 0)
        self.assertEqual(result['fn'], 1)
        self.assertNotEqual(result['smin'], float('inf'))


if __name__ == '__main__':
    unittest.main()
