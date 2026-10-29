import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]

class ConservativeCLITests(unittest.TestCase):
    def test_canonicalized_cli_export_and_audit(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/'go.obo').write_text('''format-version: 1.2

[Term]
id: GO:0032991
name: protein-containing complex
namespace: cellular_component

[Term]
id: GO:0099001
alt_id: GO:0099002
name: synthetic two-member complex
namespace: cellular_component
is_a: GO:0032991 ! protein-containing complex

''')
            (root/'input.tsv').write_text('a\tGO:0099002|0.9\nb\n')
            (root/'evidence.json').write_text(json.dumps([dict(protein='b',term='GO:0099002',source='test evidence')]))
            subprocess.run([sys.executable,str(ROOT/'pipeline/conservative_repair.py'),
                '--predictions',str(root/'input.tsv'),'--go',str(root/'go.obo'),'--threshold','0.3',
                '--promotion-evidence',str(root/'evidence.json'),'--output',str(root/'out.tsv'),
                '--audit',str(root/'audit.json')],check=True,capture_output=True,text=True)
            audit=json.loads((root/'audit.json').read_text())
            self.assertEqual(audit['normalization']['alternative_id_pairs'],1)
            self.assertEqual(audit['normalization']['hierarchy_added_pairs'],1)
            self.assertEqual(audit['demotions'],[])
            self.assertEqual({tuple(x) for x in audit['promotions']},{('GO:0099001','b'),('GO:0032991','b')})
            self.assertNotIn('GO:0099002',(root/'out.tsv').read_text())

if __name__=='__main__':
    unittest.main()
