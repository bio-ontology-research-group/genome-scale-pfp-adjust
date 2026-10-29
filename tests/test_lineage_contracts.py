"""Union ambiguity must not be reported as a uniquely inferred lineage."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.audit_lineage_results import supported

class IdentityGO:
    def get_ancestors(self, g):
        return {g}

class LineageTests(unittest.TestCase):
    def test_union_alternatives_do_not_become_forced(self):
        isa = {'NCBITaxon_2': {'NCBITaxon_1'}, 'NCBITaxon_3': {'NCBITaxon_1'}}
        only = {'g': ['NCBITaxon_Union_9']}
        unions = {'NCBITaxon_Union_9': {'NCBITaxon_2', 'NCBITaxon_3'}}
        ranks = {'NCBITaxon_2': 'domain', 'NCBITaxon_3': 'domain'}
        self.assertEqual(supported({'g'}, only, {}, isa, {}, unions, IdentityGO(), ranks), set())
        self.assertEqual(supported({'g', 'h'}, only, {'h': ['NCBITaxon_3']}, isa, {}, unions, IdentityGO(), ranks), {'NCBITaxon_2'})

if __name__ == '__main__':
    unittest.main()
