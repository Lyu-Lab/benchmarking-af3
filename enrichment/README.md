# Enrichment

## Installation
Install two conda environments `py3d` and `dock37_py27` via yaml files `./envs/py3d.yml` and `./envs/dock37_py27.yml` with
```bash
conda env create -f ./envs/dock37_py27.yml
conda env create -f ./envs/py3d.yml
```
For more information please see conda docs linked [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

The DOCK3 source code can be downloaded [here](https://github.com/docking-org/dock3) for `plots.py` and `enrich.py` and initialize `DOCKBASE` path.

```bash
export DOCKBASE=/path/to/dock37/DOCK-3.7-trunk
```

## Run analysis code
0. (optional) Set up conda `source "$CONDA_SH"` if applicable.
1. [under `py3d`] Run `split_running_sum.py` to separate the raw data file containing `metric` values for all targets into individual folders organized by target. The output follows the structure `scoring_metric_name/target_name/split.csv`. For example, running it on `DUDEZ/bootstrap_running_sum` from the data deposit linked [here](https://app.globus.org/file-manager?origin_id=67252878-c56c-45f2-9eeb-dd7f45e01720&origin_path=%2FDUDEZ%2Fbootstrap_running_sum%2F) will produce `scoring_metric_name/target_name/split.csv` files for each DUDEZ target. You can also find an example locally in `enrichment/example_running_sum` and run `python split_running_sum.py --input-dir ./example_running_sum --out-dir ./example_running_sum`.
2. Run `compute_auc_from_splits.sh` to calculate AUC and logAUC values for each DUDEZ target by setting root to the output of step 1.
For example, `bash 'compute_auc_from_splits.sh' --root ./example_running_sum/docking_score` where `docking_score` contains subfolders `target_name1`, `target_name2`...
3. [under `py3d`] Run `collect_auc_summary.py` to generate a summary `.csv` file containing AUC, logAUC, number of ligands, and number of decoys for each target.
For example,
```bash
cd ./example_running_sum
python ../collect_auc_summary.py
```
