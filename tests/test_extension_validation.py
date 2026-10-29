import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.validate_extension_exports import validate, sha

class GO:
    def get_ancestors(self,g):
        return {g,'root'}

class ExportValidationTests(unittest.TestCase):
    def run_case(self, extra='', unsupported=False):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); inp=root/'deepgo-se/cc/stage1'
            inp.mkdir(parents=True)
            source=inp/'optimized_fold_01_taxon_42.tsv'
            score=0 if unsupported else .1
            source.write_text(f'a\tc|0.9\troot|0.9\nb\tc|{score}\troot|{score}\n')
            dest=root/'conservative/deepgo-se/42.tsv';dest.parent.mkdir(parents=True)
            dest.write_text('a\tc|0.9\troot|0.9'+extra+'\nb\tc|0.301\troot|0.301\n')
            result=dest.with_suffix('.json')
            result.write_text(json.dumps(dict(predictor='deepgo-se',taxon='42',threshold=.3,
                input_sha256=sha(source),export_sha256=sha(dest),promotions=[['c','b'],['root','b']],
                demotions=[],unresolved=[],initial_singletons=1,evidence_policy='nonzero')))
            return validate(SimpleNamespace(results=root),result,GO(),{'c':{'root'}},{'c'},
                            {},{},{},{},{'b':'P12345'},{'42::c':{'P12345'}})

    def test_complete_export_is_validated(self):
        result=self.run_case()
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['membership']['listed_promotions'],1)

    def test_unrelated_score_mutation_fails(self):
        with self.assertRaises(AssertionError):
            self.run_case(extra='\tx|0.01')

    def test_unsupported_promotion_fails(self):
        with self.assertRaises(AssertionError):
            self.run_case(unsupported=True)

if __name__=='__main__':
    unittest.main()
