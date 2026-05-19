# SPROF-GO wrapper

`run_sprof_go.py` shells out to the upstream SPROF-GO inference code, then
`process_sprof-go_results.py` converts the raw outputs to the per-protein TSV
format consumed by `pipeline/run_adjustment_pipeline.py`.

## Setup

```
git clone https://github.com/biomed-AI/SPROF-GO ${SPROF_GO_REPO}
cd ${SPROF_GO_REPO}
# follow upstream README to install dependencies and download the pretrained model
```

Set `SPROF_GO_REPO=` in your top-level `config.env`.

## Reference

Yuan, Q., Tian, C., Yang, Y. (2023). Fast and accurate protein function
prediction from sequence through pretrained language model and homology-based
label diffusion. *Briefings in Bioinformatics*.
