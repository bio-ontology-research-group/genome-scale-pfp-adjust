"""Exercise complete-result gates and numerical/TeX generation with synthetic data."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('tables',ROOT/'reproducibility/reproduce_all_tables_and_figures.py')
tables=importlib.util.module_from_spec(spec); spec.loader.exec_module(tables)


def fixture(base):
    (base/'summaries').mkdir()
    inputs={}; checks=[]; summaries=[]
    for m in tables.MODELS:
        n=219 if m.endswith('_heldout') else 50
        for a in tables.ASPECTS:
            sources=[]
            for t in map(str,range(n)):
                h=f'{m}/{a}/{t}'
                sources.append({'task':[m,a,t],'result_sha256':h,'code':{'runner':'fixed'},'configuration':{'solver_seconds':120}})
                inputs[h]=h
                if m in tables.MODELS[:3] and a=='cc':
                    for mode in ['standard','demotion']:
                        checks.append({'predictor':m,'taxon':t,'mode':mode,'metadata_sha256':h,
                            'stage1_sha256':'input','stage2_sha256':mode,'promotions':1,'demotions':2})
            native=.3 if m=='deepgo-se' else .001
            conditions={(mode,t) for mode in ['before','stage1'] for t in ([native] if m=='seq_sim_heldout' else [.001,.01,.05,.1,.2,.3,.4,.5])}
            conditions.add(('stage1_unknown',native))
            if m in tables.MODELS[:3]:
                conditions.add(('published_funtaxis',native))
                if a=='cc':
                    conditions.update((mode,native) for mode in ['stage2_standard','stage2_demotion']+[f'soft_stage{s}_{p}' for p in [.2,.5,1.0] for s in [1,2]])
            rows=[]
            for mode,t in sorted(conditions):
                rows.append({'predictor':m,'aspect':a,'mode':mode,'tau':t,'n_organisms':n,
                    'n_annotated_proteins':n,'n_annotated_organisms':n,'active_input':1000,'promotions':n,'demotions':2*n,
                    'demotion_active_pct':n/5,'flip_grid_pct':.1,'known_taxon_consistency':1,'gaef_taxon_consistency':1,
                    'scope_counts':{'n_mask_proteins':n,'n_annotation_matches':n},
                    'identifier_normalization':{},
                    'metrics':{k:.5 for k in ['fmax','smin','avg_auc','aupr','precision','recall']},
                    'complex_coherence':1 if a=='cc' else None,'process_coherence':1,'process_coherence_n':n,
                    'pathway_coherence':None,'pathway_coherence_n':0,'essential_completeness':None,'essential_completeness_n':0,
                    'solver_status_counts':{'OPTIMAL':n},'solver_objective_sum':1,'solver_bound_sum':1})
            (base/'summaries'/f'{m}_{a}.json').write_text(json.dumps({'schema_version':2,'sources':sources,'rows':rows,
                'ic_reference_inputs':{'annotations':'fixed'},'ic_split_inputs':{'split':'fixed'},
                'reference_manifest_sha256':'fixed-reference'}))
        for a in tables.ASPECTS+['joint']:
            summaries.append({'predictor':m,'aspect':a,'n':n,'supported':{'deepest_correct_rank':{'domain':n},'first_divergence_rank':{'none':n},'exact_prefix_count':n},
                'cross_aspect_rank_agreement':{r:{'compatible':n,'all_three_resolved_equal':0} for r in tables.RANKS[1:]}})
    (base/'lineage.json').write_text(json.dumps({'summary':summaries,'input_sha256':inputs}))
    (base/'export-validation.json').write_text(json.dumps({'validated':300,'records':checks}))
    (base/'evidence-sensitivity.json').write_text(json.dumps({'rows':[{'taxon':t,'complex_counts':{'curated':{'present':2,'singletons':1},'all_proteome':{'present':3,'singletons':0}},'curated_singletons_with_partner_all_same_proteins':1,'curated_singletons_with_partner_all_proteome':1} for t in range(50)]}))
    for m in tables.MODELS[:3]:
        for mode in ['standard','demotion']:
            (base/f'membership_{m}_{mode}.json').write_text(json.dumps({'n_organism_pairs':50,
                'provenance':{'pairs':[{'taxon':str(t),'stage1':'input','stage2':mode} for t in range(50)]},
                'promotions':{'total':50,'with_curated_reference':2,'true_subunit_promotions':1},
                'demotions':{'total':100,'with_curated_reference':2,'curated_subunit_demotions':1},
                'unmapped_promotions_with_reference':0,'unmapped_demotions_with_reference':0}))


class TableTests(unittest.TestCase):
    def test_empty_denominators_render_as_undefined(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as d:
            base=Path(d); fixture(base)
            path=base/'summaries/seq_sim_timeset_cc.json'
            data=json.loads(path.read_text())
            row=next(r for r in data['rows'] if r['mode']=='stage1' and r['tau']==.001)
            row.update(active_input=0,demotions=0,demotion_active_pct=None,flip_grid_pct=None)
            path.write_text(json.dumps(data))
            tables.generate(argparse.Namespace(results=base,output=base/'generated'))
            with (base/'generated/main-flips.tsv').open() as fh:
                rows=list(csv.reader(fh,delimiter='\t'))
            self.assertEqual(rows[1][-2:],['--','--'])

    def test_complete_bundle_and_rejection_of_stale_validation(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as d:
            base=Path(d); fixture(base); out=base/'generated'
            args=argparse.Namespace(results=base,output=out)
            tables.generate(args)
            text=(out/'main-results.tex').read_text()
            self.assertIn(r'$F_{\max}$',text)
            self.assertIn(r'\ref{tab:main-metrics}',text)
            for path in out.glob('*.tex'):
                self.assertNotIn(b'\r',path.read_bytes())
            with (out/'supp-threshold-structure.tsv').open() as fh:
                self.assertEqual(len(list(csv.reader(fh,delimiter='\t'))),97)
            manifest=out/'generation-manifest.json'; old=manifest.read_bytes()
            path=base/'export-validation.json'; data=json.loads(path.read_text()); data['records'][0]['metadata_sha256']='stale'; path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'different result metadata'):
                tables.generate(args)
            self.assertEqual(manifest.read_bytes(),old)
            path.unlink()
            with self.assertRaises(FileNotFoundError):
                tables.generate(args)
            self.assertEqual(manifest.read_bytes(),old)

if __name__=='__main__': unittest.main()
