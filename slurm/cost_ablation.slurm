#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=cost_ablation
#SBATCH --output=%x.log
#
# Cost-function ablation: re-run Stage-1 (taxon-guided) on DeepGO-SE CC under
# three cost forms (margin = above-threshold mass [default], uniform, IC) and
# report the flip rate and pairwise flip-set agreement, to show the qualitative
# result is unchanged by the cost choice.
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

ANN=${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic
BASE=${DATA_DIR}/swissprot_proteomes_folds/deepgo-se_results/cc/predictions
OUTROOT=${OUT_DIR}/cost_ablation
IC=${OUTROOT}/go_ic.tsv
RES=${OUTROOT}/fliprates.tsv
mkdir -p ${OUTROOT}; : > ${RES}

[ -f ${IC} ] || python analysis/compute_go_ic.py --annotations_dir ${ANN} --output_file ${IC}

export PFP_IC_FILE=${IC}
for MODE in margin uniform ic; do
    export PFP_COST_MODE=${MODE}
    OD=${OUTROOT}/${MODE}
    echo "=== cost mode ${MODE} ==="
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE} --optimized_dir ${OD} \
        --optimized --default_threshold 0.3 --provide_taxon_id 2>&1 | tail -2
    python analysis/flip_rate.py --predictions_dir ${BASE} --optimized_dir ${OD} \
        --threshold 0.3 --label deepgo-se_cc_${MODE} | tee -a ${RES}
done

echo "=== pairwise flip-set agreement (Jaccard over demoted (protein,term)) ==="
python3 - <<PYEOF
import glob, os, re, itertools, json
def demoted(pred_dir, opt_dir, tau=0.3):
    fn=re.compile(r'_fold_(\d+)_taxon_(\d+)\.tsv$'); S=set()
    def above(p):
        d={}
        for l in open(p):
            x=l.rstrip('\n').split('\t')
            if not x or not x[0]: continue
            d[x[0]]=set(t.split('|')[0] for t in x[1:] if '|' in t and float(t.split('|')[1])>tau)
        return d
    for opt in glob.glob(os.path.join(opt_dir,'optimized_fold_*_taxon_*.tsv')):
        m=fn.search(os.path.basename(opt));
        if not m: continue
        pr=os.path.join('${BASE}','predictions_fold_%s_taxon_%s.tsv'%(m.group(1),m.group(2)))
        if not os.path.exists(pr): continue
        a,b=above(pr),above(opt)
        for prot,ts in a.items():
            for go in ts-b.get(prot,set()): S.add((m.group(2),prot,go))
    return S
modes=['margin','uniform','ic']
sets={m:demoted('${BASE}', '${OUTROOT}/'+m) for m in modes}
for a,b in itertools.combinations(modes,2):
    A,B=sets[a],sets[b]; j=len(A&B)/len(A|B) if (A|B) else 1.0
    print('%s vs %s: |A|=%d |B|=%d Jaccard=%.4f'%(a,b,len(A),len(B),j))
PYEOF
echo "DONE cost_ablation"
