# Prospective

## Run Boltz

You can run Boltz via command `boltz predict input_path` where input_path is `boltz_input_example.yaml`.

See further Boltz documentation [here](https://github.com/jwohlwend/boltz)

## Correlation
1. Create a minimal conda environment using `pip install -r requirements.txt` or existing environment with specified packages.
2. Download `cpm_ic50.csv` which is stored in the separate data deposit [here](https://app.globus.org/file-manager?origin_id=db5285b8-e50a-4b07-94eb-890c1e27d127&origin_path=%2Fsigma2_prospective%2F) and run the Spearman correlation between the experimental log IC50 value and predicted affinity value and affinity probability.

```bash
cd prospective

python correlation.py --experimental_data_path cpm_ic50.csv --boltz2_data_path boltz2_data.csv
```

Installation and running of the demo should be under 10 minutes.