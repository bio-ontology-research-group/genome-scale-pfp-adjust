"""
Local patches/extensions to the upstream GAEF package.

Upstream GAEF (https://github.com/bio-ontology-research-group/GAEF) is used
unmodified for `GAEF.completeness`, `GAEF.utils.Ontology`, `FUNC_DICT`, and
`NAMESPACES`. Everything that diverges from upstream lives here:

  - coherence              : adds `include_part_of` kwarg to parse_go_ontology,
                             switches has_part to genome-level + new return shape.
  - taxon_consistency      : Python port of taxon_consistency.groovy.
  - complex_classifier     : Python port of complex coherence helpers.
  - Ontology IC methods    : `calculate_ic` / `get_ic` / `get_norm_ic` are
                             attached to the upstream Ontology class on import,
                             matching what the local GAEF fork carries.

Vendored constraint files live in <repo>/data/constraints/.

If/when these patches are upstreamed, this package can be shrunk or removed
and imports can revert to `from GAEF.coherence import ...` etc.
"""

from collections import Counter
import math

from GAEF.utils import Ontology as _UpstreamOntology


def _calculate_ic(self, annotations):
    """Resnik IC: -log2(p(term)), p = count/total. Stores in self.ic; sets self.ic_norm."""
    cnt = Counter()
    total_annotations = 0
    for annot_set in annotations:
        for term in annot_set:
            if term in self.ont:
                cnt[term] += 1
                total_annotations += 1
    self.ic = {}
    self.ic_norm = 0.0
    if total_annotations == 0:
        return
    for go_id, count in cnt.items():
        ic_value = -math.log2(count / total_annotations)
        self.ic[go_id] = ic_value
        self.ic_norm = max(self.ic_norm, ic_value)


def _get_ic(self, go_id):
    if self.ic is None:
        raise Exception("IC not yet calculated. Call calculate_ic() first.")
    return self.ic.get(go_id, 0.0)


def _get_norm_ic(self, go_id):
    if self.ic_norm == 0.0:
        return 0.0
    return self.get_ic(go_id) / self.ic_norm


# Attach idempotently so re-imports don't double-bind.
if not hasattr(_UpstreamOntology, "calculate_ic"):
    _UpstreamOntology.calculate_ic = _calculate_ic
if not hasattr(_UpstreamOntology, "get_ic"):
    _UpstreamOntology.get_ic = _get_ic
if not hasattr(_UpstreamOntology, "get_norm_ic"):
    _UpstreamOntology.get_norm_ic = _get_norm_ic
