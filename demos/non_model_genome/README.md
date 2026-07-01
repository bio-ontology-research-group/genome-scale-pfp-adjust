# Non-model genome demo: Thermococcus kodakarensis

*2026-07-01T05:08:00Z by Showboat 0.6.1*
<!-- showboat-id: c88e60e8-b90e-4967-83c8-42c7def302b6 -->

This executable demo runs GenoAdjust on a small non-model archaeal genome example: *Thermococcus kodakarensis* KOD1 (`NCBITaxon_69014`). The input is intentionally tiny so it can be run on a laptop, but it uses the same `pipeline/run_adjustment_pipeline.py`, Stage 1 taxon-consistency solver, and Stage 2 complex-coherence solver used by the paper workflow.

The demo prediction file contains six protein-like records. Two annotations are biologically incompatible with an archaeon lineage: `GO:0045087` (*innate immune response*) and `GO:0032501` (*multicellular organismal process*), both marked as `never_in_taxon NCBITaxon_2157` in the miniature constraint file. The same input also contains `GO:0009333` (*cysteine synthase complex*) above threshold on exactly one protein, so Stage 2 must repair the singleton complex by either demoting it or promoting a partner.

Run all commands from the repository root. The generated files under `demos/non_model_genome/output/` are ignored by git and can be regenerated at any time.

First confirm the lightweight solver dependencies are available. The command uses `uv run` to provide `ortools` and `tqdm` in an isolated environment; users who already installed `requirements-test.txt` can run the Python scripts directly.

```bash
uv run --quiet --with ortools==9.15.6755 --with tqdm python - <<'PY'
from ortools.sat.python import cp_model
import tqdm
print("Dependencies available: ortools + tqdm")
PY

```

```output
Dependencies available: ortools + tqdm
```

The prediction input follows the repository's standard per-genome TSV format: first column is the protein identifier; each remaining column is `GO_ID|score`. The filename includes the fold and species taxon so the pipeline can discover it automatically.

```bash
sed -n '1,20p' demos/non_model_genome/data/predictions/predictions_fold_00_taxon_69014.tsv

```

```output
TK0001	GO:0045087|0.92	GO:0006412|0.77	GO:0009333|0.86
TK0002	GO:0032501|0.88	GO:0005737|0.70	GO:0009333|0.28
TK0003	GO:0006412|0.74	GO:0009333|0.12
TK0004	GO:0005737|0.68
TK0005	GO:0002181|0.55	GO:0006412|0.52
TK0006	GO:0000166|0.62	GO:0003677|0.61
```

The miniature constraint set is deliberately scoped to this demo. It maps *T. kodakarensis* to Archaea (`NCBITaxon_69014 -> NCBITaxon_2157`) and includes only the GO/taxon and complex hierarchy edges needed to exercise the two solver stages offline.

```bash
python3 - <<'PY'
from pathlib import Path
base = Path('demos/non_model_genome/data/constraints')
for name in ['go_taxon_constraints.tsv', 'taxon_hierarchy.tsv', 'protein_complexes.tsv']:
    print(f'## {name}')
    print((base / name).read_text().strip())
PY

```

```output
## go_taxon_constraints.tsv
GO_ID	Constraint_Type	Taxon_ID
GO:0032501	never_in_taxon	NCBITaxon_2157
GO:0045087	never_in_taxon	NCBITaxon_2157
GO:0006412	only_in_taxon	NCBITaxon_131567
GO:0005737	only_in_taxon	NCBITaxon_131567
## taxon_hierarchy.tsv
Term	Relationship	Parent/Disjoint_From_Term
NCBITaxon_69014	is_a	NCBITaxon_2157
NCBITaxon_2157	is_a	NCBITaxon_131567
## protein_complexes.tsv
GO_term	classification	definition
GO:0009333	n	Cysteine synthase is a multienzyme complex made up of serine acetyltransferase and O-acetylserine (thiol)-lyase subunits.
```

Run the full two-stage adjustment pipeline. Stage 1 receives the known organism taxon (`--provide_taxon_id`), and Stage 2 uses the optimized complex-coherence formulation with `top_k=3`, enough to include the existing singleton and the two closest candidate partners in this small example.

```bash
set -euo pipefail
rm -rf demos/non_model_genome/output
uv run --quiet --with ortools==9.15.6755 --with tqdm python pipeline/run_adjustment_pipeline.py \
  --predictions_dir demos/non_model_genome/data/predictions \
  --optimized_dir demos/non_model_genome/output \
  --fold_thresholds demos/non_model_genome/data/fold_thresholds.json \
  --constraints_file demos/non_model_genome/data/constraints/go_taxon_constraints.tsv \
  --go_hierarchy_file demos/non_model_genome/data/constraints/go_hierarchy.tsv \
  --taxon_hierarchy_file demos/non_model_genome/data/constraints/taxon_hierarchy.tsv \
  --ncbitaxon_hierarchy_file demos/non_model_genome/data/constraints/ncbitaxon_hierarchy.tsv \
  --complexes_file demos/non_model_genome/data/constraints/protein_complexes.tsv \
  --provide_taxon_id \
  --complex_coherence \
  --optimized \
  --top_k 3 > demos/non_model_genome/output.log 2>&1
python3 demos/non_model_genome/check_demo_output.py

```

```output
Demo genome: Thermococcus kodakarensis KOD1 (NCBITaxon_69014)
Proteins: 6
Above-threshold annotations before: 11
Above-threshold annotations after: 10
Taxon-stage removals:
  - GO:0032501 (multicellular organismal process) from TK0002
  - GO:0045087 (innate immune response) from TK0001
Complex term GO:0009333 before: TK0001
Complex term GO:0009333 after: TK0001, TK0002
OK: taxon consistency and complex coherence invariants hold.
```

The adjusted TSV shows the concrete edits. `GO:0045087` and `GO:0032501` are gone. The singleton `GO:0009333` remains on `TK0001`, and `TK0002` is promoted from `0.28` to `0.301`, just above the `0.3` threshold, because that is the minimum-cost complex-coherence repair.

```bash
python3 - <<'PY'
from pathlib import Path
path = Path('demos/non_model_genome/output/optimized_fold_00_taxon_69014.tsv')
for line in path.read_text().splitlines():
    protein, *terms = line.split('\t')
    print('\t'.join([protein] + sorted(terms)))
PY

```

```output
TK0001	GO:0006412|0.770000	GO:0009333|0.860000
TK0002	GO:0005737|0.700000	GO:0009333|0.301000
TK0003	GO:0006412|0.740000	GO:0009333|0.120000
TK0004	GO:0005737|0.680000
TK0005	GO:0002181|0.550000	GO:0006412|0.520000
TK0006	GO:0000166|0.620000	GO:0003677|0.610000
```

To re-run this proof document, execute `uvx showboat --workdir "$PWD" verify demos/non_model_genome/README.md` from the repository root. Verification reruns every command above and diffs the captured output.

