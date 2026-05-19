# DeepGO-SE wrapper

`predict_deepgo-se.py` is a thin wrapper around the upstream DeepGO-SE model.
The model weights and inference code live in a separate repo.

## Setup

```
# clone next to genome-scale-pfp-adjust
git clone https://github.com/bio-ontology-research-group/deepgo2 ${DEEPGO_SE_REPO}
cd ${DEEPGO_SE_REPO}
# download the pretrained checkpoints
wget https://deepgo.cbrc.kaust.edu.sa/data/deepgo2/data.tar.gz
tar -xzf data.tar.gz
```

Set `DEEPGO_SE_REPO=` in your top-level `config.env` to point at the clone, and
make sure its checkpoint dir is reachable from this script's CLI args.

## Reference

Kulmanov, M., Hoehndorf, R. (2022). DeepGO-SE: Protein function prediction as
approximate semantic entailment. *Nature Machine Intelligence*.
