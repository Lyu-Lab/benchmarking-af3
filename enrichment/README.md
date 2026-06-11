# Enrichment

## Installation
Install two conda environments, `py3d` and `dock37_py27`.

**Linux** — use the pinned cluster exports:
```bash
conda env create -f ./envs/dock37_py27.yml
conda env create -f ./envs/py3d.yml
```

**macOS (Apple Silicon or Intel)** — use the minimal env files:
```bash
bash ./setup_mac_envs.sh
```
Or manually:
```bash
conda env create -f ./envs/py3d_osx.yml
CONDA_SUBDIR=osx-64 conda env create -f ./envs/dock37_py27_osx.yml   # Apple Silicon only; omit on Intel Mac
```
On Apple Silicon, `dock37_py27` runs under Rosetta because Python 2.7 is not available for `osx-arm64`.

For more information please see conda docs linked [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)


## Run analysis code
0. (optional) Set up conda `source "$CONDA_SH"` if applicable.
1. [under `py3d`] Run `split_running_sum.py` to separate the raw data file containing `metric` values for all targets into individual folders organized by target. The output follows the structure `scoring_metric_name/target_name/split.csv`. For example, running it on `DUDEZ/bootstrap_running_sum` from the data deposit linked [here](https://app.globus.org/file-manager?origin_id=db5285b8-e50a-4b07-94eb-890c1e27d127&origin_path=%2FDUDEZ%2Fbootstrap_running_sum%2F) will produce `scoring_metric_name/target_name/split.csv` files for each DUDEZ target. You can also find an example locally in `enrichment/example_running_sum` and run `python split_running_sum.py --input-dir ./example_running_sum --out-dir ./example_running_sum`.
2. Run `compute_auc_from_splits.sh` to calculate AUC and logAUC values for each DUDEZ target by setting root to the output of step 1.
For example, `bash 'compute_auc_from_splits.sh' --root ./example_running_sum/docking_score` where `docking_score` contains subfolders `target_name1`, `target_name2`...
3. [under `py3d`] Run `collect_auc_summary.py` to generate a summary `.csv` file containing AUC, logAUC, number of ligands, and number of decoys for each target.
For example,
```bash
cd ./example_running_sum
python ../collect_auc_summary.py
```

Installation and running of the demo should be under 10 minutes.