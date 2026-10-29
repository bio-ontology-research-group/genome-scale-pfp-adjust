"""Source labels replace inherited labels without changing split identity."""
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from data_prep.prepare_baseline_inputs import prepare, digest

class BaselineInputTests(unittest.TestCase):
    def test_source_labels_replace_unsafe_ancestors_and_empty_rows(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as root:
            root=Path(root); refs=root/'refs'; refs.mkdir(); terms=root/'terms'; terms.mkdir()
            text='format-version: 1.2\n\n'
            for aspect,ns,term in [('cc','cellular_component','GO:0005575'),('mf','molecular_function','GO:0003674'),('bp','biological_process','GO:0008150')]:
                text+=f'[Term]\nid: {term}\nname: {aspect}\nnamespace: {ns}\n\n'
                (terms/aspect).mkdir();pd.DataFrame({'gos':[term]}).to_pickle(terms/aspect/'terms.pkl')
            go=root/'go.obo';go.write_text(text)
            ann=refs/'annots_taxon_1.tsv';ann.write_text('p\tGO:0005575\nq\tGO:0003674\nr\t\n')
            (refs/'reference-manifest.json').write_text(json.dumps({'go_sha256':digest(go),'records':[{'taxon':'1','reference_file':ann.name,'reference_sha256':digest(ann)}]}))
            train=root/'train.tsv';train.write_text('p\t1\tGO:0005575\tGO:0008150\nr\t1\tGO:0008150\n')
            val=root/'val.tsv';val.write_text('q\t1\tGO:0003674\n')
            test=root/'test.txt';test.write_text('2\n')
            args=SimpleNamespace(references=refs,go=go,train=train,validation=val,test_taxa=test,terms_root=terms,output=root/'out')
            result=prepare(args)
            self.assertEqual(result['split_counts']['train'],dict(original=2,retained=1,empty_removed=1,changed=2))
            self.assertEqual((root/'out/train_proteins.tsv').read_text(),'p\t1\tGO:0005575\n')
            args.output=root/'bad';val.write_text('p\t1\tGO:0005575\n')
            with self.assertRaisesRegex(ValueError,'overlap'):
                prepare(args)

if __name__=='__main__': unittest.main()
