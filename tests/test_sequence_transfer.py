"""Sequence-transfer source and score contracts, without running an aligner."""
import gzip
import json
from types import SimpleNamespace
import pandas as pd
import os
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from baselines.seq_sim.rebuild_from_sources import sequences, transfer_scores, run
from data_prep.prepare_baseline_inputs import digest

class SequenceTests(unittest.TestCase):
    def test_bitscore_weighting_uses_all_annotated_hits(self):
        scores=transfer_scores({'a':3,'b':1},{'a':{'parent','child'},'b':{'parent','other'}})
        self.assertEqual(scores,{'parent':1,'child':.75,'other':.25})
        self.assertEqual(transfer_scores({},{}),{})
        with self.assertRaises(ValueError):transfer_scores({'a':float('nan')},{'a':{'term'}})
        with self.assertRaises(ValueError):transfer_scores({'unknown':1},{})

    def test_sequences_require_complete_entries(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as root:
            path=Path(root)/'source.dat.gz'
            with gzip.open(path,'wt') as handle:
                handle.write('ID   A_X Reviewed; 4 AA.\nSQ   SEQUENCE 4 AA;\n     MA AA\n//\n')
            self.assertEqual(list(sequences(path)),[('A_X','MAAA')])
            with gzip.open(path,'wt') as handle:
                handle.write('ID   A_X Reviewed; 4 AA.\nSQ   SEQUENCE 4 AA;\n     MAAA\n')
            with self.assertRaisesRegex(ValueError,'Truncated'):list(sequences(path))

    def test_complete_query_universe_includes_proteins_without_hits(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            root=Path(folder);refs=root/'refs';refs.mkdir();terms=root/'terms';terms.mkdir()
            source=root/'source';source.mkdir();records=[]
            go=root/'go.obo';text='format-version: 1.2\n\n'
            for aspect,ns,term in [('cc','cellular_component','GO:0005575'),('mf','molecular_function','GO:0003674'),('bp','biological_process','GO:0008150')]:
                text+=f'[Term]\nid: {term}\nname: {aspect}\nnamespace: {ns}\n\n'
                (terms/aspect).mkdir();pd.DataFrame({'gos':[term]}).to_pickle(terms/aspect/'terms.pkl')
            go.write_text(text)
            for taxon,proteins in [('1',['A']),('2',['B','C'])]:
                path=source/'Bacteria'/taxon/(taxon+'.dat.gz');path.parent.mkdir(parents=True)
                with gzip.open(path,'wt') as handle:
                    for protein in proteins:handle.write(f'ID   {protein} Reviewed; 1 AA.\nSQ   SEQUENCE 1 AA;\n     M\n//\n')
                ref=refs/f'annots_taxon_{taxon}.tsv'
                ref.write_text(''.join(p+'\tGO:0005575\n' for p in proteins))
                records.append(dict(taxon=taxon,source_file=str(path),source_sha256=digest(path),reference_file=ref.name,reference_sha256=digest(ref)))
            (refs/'reference-manifest.json').write_text(json.dumps(dict(go_sha256=digest(go),records=records)))
            train=root/'train.txt';train.write_text('1\n');test=root/'test.txt';test.write_text('2\n')
            diamond=root/'diamond'
            diamond.write_text('#!'+sys.executable+'\nimport sys\nfrom pathlib import Path\na=sys.argv[1:]\nif a[0]=="version": print("fixture aligner")\nelif a[0]=="makedb": Path(a[a.index("-d")+1]+".dmnd").write_text("fixture")\nelse: Path(a[a.index("-o")+1]).write_text("B\\tA\\t100\\t1e-20\\n")\n')
            diamond.chmod(0o755)
            args=SimpleNamespace(output=root/'out',references=refs,go=go,train_taxa=train,test_taxa=test,terms_root=terms,source_root=source,diamond=str(diamond),threads=1,score_floor=.001)
            run(args)
            rows=(root/'out/cc/predictions/predictions_fold_01_taxon_2.tsv').read_text().splitlines()
            self.assertEqual(rows,['B\tGO:0005575|1','C'])
            manifest=json.loads((root/'out/prediction-manifest.json').read_text())
            self.assertEqual(manifest['query_proteins'],2)
            self.assertEqual(manifest['queries_with_hits'],1)

if __name__=='__main__':unittest.main()
