# Enrichment

## Installation
Install conda environment `py3d`.

**Linux** — use the pinned cluster exports:
```bash
conda env create -f ./envs/py3d.yml
```

**macOS (Apple Silicon or Intel)** — use the minimal env files:
```bash
conda env create -f ./envs/py3d_osx.yml
```

For more information please see conda docs linked [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

## Run analysis code
0. (optional) Set up conda `source "$CONDA_SH"` if applicable.
1. Run `split_running_sum.py` to separate the raw data file containing `metric` values for all targets into individual folders organized by target. The output follows the structure `scoring_metric_name/target_name/split.csv`. For example, running it on `DUDEZ/bootstrap_running_sum` from the data deposit linked [here](https://app.globus.org/file-manager?origin_id=db5285b8-e50a-4b07-94eb-890c1e27d127&origin_path=%2FDUDEZ%2Fbootstrap_running_sum%2F) will produce `scoring_metric_name/target_name/split.csv` files for each DUDEZ target. You can also find an example locally in `enrichment/example_running_sum` and run `python split_running_sum.py --input-dir ./example_running_sum --out-dir ./example_running_sum`.
2. Run `tie_safe_auc.py` to calculate AUC and logAUC values for each target by setting root to the output of step 1.
For example, `python3 tie_safe_auc.py example_running_sum/docking_score/AMPC/split.csv --scores-col-name docking_score` would produce the output AUC=76.36 logAUC=14.59 in bash.

Installation and running of the demo should be under 10 minutes.
