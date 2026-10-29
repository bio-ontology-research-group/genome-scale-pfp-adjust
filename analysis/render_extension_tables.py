"""Generate supplementary extension tables from the checked summary only."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.revision_extensions import sha

MODELS={'seq_sim_timeset':'SS-E','deepgo-se':'DGO-E','sprof-go':'SP-E',
        'seq_sim_heldout':'SS-H','mlp_heldout':'MLP-H'}
ORDER=list(MODELS)


def table(path,label,caption,columns,rows):
    lines=[r'\begingroup\scriptsize',r'\setlength{\tabcolsep}{4pt}',
           r'\begin{longtable}{'+''.join('l' if i<2 else 'r' for i in range(len(columns)))+'}',
           r'\caption{'+caption+r'}\label{tab:'+label+r'}\\',r'\toprule',
           ' & '.join(columns)+r' \\',r'\midrule',r'\endfirsthead',r'\toprule',
           ' & '.join(columns)+r' \\',r'\midrule',r'\endhead']
    lines.extend(' & '.join(map(str,row))+r' \\' for row in rows)
    lines.extend([r'\bottomrule',r'\end{longtable}',r'\endgroup'])
    path.write_text('\n'.join(lines)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();data=json.loads(args.summary.read_text())
    assert data['counts']==dict(lineage_tasks=588,lineage_outcomes=2352,rank_rows=16464,validated_exports=300,recovery_tasks=150)
    args.output.mkdir(parents=True,exist_ok=True)
    rows=sorted(data['lineage'],key=lambda x:(ORDER.index(x['predictor']),['cc','mf','bp','joint'].index(x['aspect'])))
    table(args.output/'supp-optimal-lineage.tex','supp-optimal-lineage',
          'Taxonomic entailment across all optimal GO repairs. Required taxa count ranked assertions across organisms. '
          'A correct prefix contains at least one required ranked taxon and agrees with the true ancestral lineage. '
          'Unknown counts undecided feasibility checks. SS: Seq-Sim; DGO: DeepGO-SE; SP: SPROF-GO; MLP: MLP-ESM2; '
          'E: entry-date split; H: heldout split.',
          ['Model','Aspect','$N$','Selected repair','All optima','Correct prefixes','Unknown'],
          [[MODELS[r['predictor']],r['aspect'].upper(),r['n'],r['conditional_claims'],r['robust_claims'],
            r['robust_correct_prefix'],r['undetermined']] for r in rows])
    modes=['hard_standard','hard_demotion','conservative','corroborated']
    names={'hard_standard':'Hard standard','hard_demotion':'Demotion-only',
           'conservative':'Nonzero gate','corroborated':'Corroboration gate'}
    rows=sorted(data['complexes'],key=lambda x:(ORDER.index(x['predictor']),modes.index(x['mode'])))
    table(args.output/'supp-conservative-complexes.tex','supp-conservative-complexes',
          'Complex repair across 50 entry-date organisms per predictor. Additions and deletions count all GO pairs. '
          'Residual gives unresolved complexes divided by initial singleton complexes. Listed/covered counts added '
          'pairs in curated membership lists among additions with an assessable reference. Listed deletions count '
          'removed curated members. Opt. gives proven optimal tasks out of 50, within each mode\'s candidate domain. '
          'Both conservative gates preserve all active inputs. Corroboration requires both other predictors above '
          'their native thresholds; nonzero uses positive primary scores. Missing membership remains unassessed.',
          ['Model','Mode','Added','Deleted','Residual','Listed/covered','Listed deletions','Opt.'],
          [[MODELS[r['predictor']],names[r['mode']],r['promotions'],r['demotions'],
            f"{r['unresolved']}/{r['initial_singletons']}",f"{r['listed']}/{r['assessable']}",
            r['listed_demotions'],f"{r['optimal']}/50"] for r in rows])
    modes=['standard','demotion_only','conservative','highest_score','conservative_corroborated','highest_score_corroborated']
    names={'standard':'Hard standard','demotion_only':'Demotion-only','conservative':'Primary-score gate',
           'highest_score':'Primary-score greedy','conservative_corroborated':'Corroboration gate',
           'highest_score_corroborated':'Corroboration greedy'}
    rows=sorted(data['recovery'],key=lambda x:(ORDER.index(x['predictor']),modes.index(x['mode'])))
    table(args.output/'supp-targeted-recovery.tex','supp-targeted-recovery',
          'Controlled completion of two-member complex annotations. Each case supplies one curated member and '
          'withholds the other from repair inputs. Other target scores are compressed below the activation '
          'threshold. Primary-score eligibility uses original scores above the native threshold; corroboration '
          'uses both other predictors. Greedy adds the eligible protein with the highest original primary score. '
          'Recovery gives hidden members restored per case. Deleted counts supplied members removed; unresolved '
          'counts singleton outputs. Other proposals count proteins apart from the supplied and hidden members. Both directions and related '
          'complex terms are correlated cases.',
          ['Model','Mode','Recovered/cases','Deleted','Unresolved','Other proposals'],
          [[MODELS[r['predictor']],names[r['mode']],f"{r['recovered']}/{r['cases']}",
            r['seed_removed'],r['unresolved'],r['other_proposals']] for r in rows])
    manifest=dict(summary_sha256=sha(args.summary),generator_sha256=sha(Path(__file__)),
                  outputs={name:sha(args.output/name) for name in ['supp-optimal-lineage.tex',
                      'supp-conservative-complexes.tex','supp-targeted-recovery.tex']})
    (args.output/'extension-generation-manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    main()
