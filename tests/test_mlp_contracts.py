"""Focused baseline regressions; requires the recorded PyTorch environment."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import torch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from baselines.mlp_esm2 import train_mlp_heldout as mlp

class MLPTests(unittest.TestCase):
    def test_features_preserve_order_and_reject_missing_or_wrong_dimension(self):
        frame=pd.DataFrame({'proteins':['a','b','c'],'orgs':['1','2','1']})
        values={'1':{'a':np.array([1.,2.]),'c':np.array([5.,6.])},'2':{'b':np.array([3.,4.])}}
        with patch.object(mlp,'load_esm2_embeddings_taxon',side_effect=lambda _,t:values[t]) as loader:
            result,missing=mlp.load_features_uniprot(frame,'unused',2)
            np.testing.assert_array_equal(result.numpy(),[[1,2],[3,4],[5,6]])
            self.assertEqual(missing,[]);self.assertEqual(loader.call_count,2)
            values['1'].pop('c')
            with self.assertRaisesRegex(ValueError,'Missing embeddings'):mlp.load_features_uniprot(frame,'unused',2)
            values['1']['c']=np.array([1.])
            with self.assertRaisesRegex(ValueError,'dimension'):mlp.load_features_uniprot(frame,'unused',2)

    def test_legacy_test_evaluation_fails_before_training(self):
        from click.testing import CliRunner
        for option in [['--no-save-predictions'], ['--annotations-go-basic-dir', 'unused']]:
            result=CliRunner().invoke(mlp.main,option)
            self.assertEqual(result.exit_code,2,result.output)
            self.assertIn('retired',result.output)
        with self.assertRaisesRegex(RuntimeError,'retired'):
            mlp.evaluate_on_annotated_test_proteins(None,[],None,[],None,None,0,None)

    def test_empty_scores_do_not_predict_at_zero_threshold(self):
        self.assertEqual(mlp.compute_fmax(np.array([[1.,0.]]),np.array([[0.,0.]]))[0],0)
        self.assertEqual(mlp.compute_fmax(np.array([[1.,0.]]),np.array([[1.,0.]]))[0],1)

    def test_train_and_export_full_precision(self):
        frame=pd.DataFrame({'proteins':['a','b'],'orgs':['1','1'],'prop_annotations':[['GO:1'],['GO:2']]})
        with patch.object(mlp,'load_features_uniprot',return_value=(torch.tensor([[1.,0.],[0.,1.]]),[])):
            model,info=mlp.train_mlp(frame,['GO:1','GO:2'],'cpu','unused',2,epochs=1,batch_size=2)
        self.assertEqual(info['n_missing_features'],0)
        self.assertTrue(all(torch.isfinite(p).all() for p in model.parameters()))
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR',str(Path.home()/'tmp'))) as folder:
            target=Path(folder)/'scores.tsv'
            score=.001000000000000001
            mlp._process_chunk_and_write((0,[('a',{'GO:1':score})],{'GO:1':{'GO:0'}},.001,str(target)))
            values=dict(t.split('|') for t in target.read_text().strip().split('\t')[1:])
            self.assertEqual(float(values['GO:1']),score)
            self.assertEqual(float(values['GO:0']),score)

if __name__=='__main__': unittest.main()
