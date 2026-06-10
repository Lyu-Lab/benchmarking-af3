# Benchmarking AlphaFold3 for structure-based ligand discovery

This repository contains the code used in our paper, "AlphaFold3 for Structure-guided Ligand Discovery."

## Overview

This repository provides pipelines, evaluation scripts, and Jupyter notebooks to run the benchmarking experiments, compare AlphaFold3 predictions with experimental structures for ligand-discovery tasks, and regenerate the manuscript's figures and tables. See the notebooks and the rest of this README for setup and reproduction instructions.

## DUDEz and LTS
The `DUDEz` and `LTS` folers contain example templates to run Alphafold3.

## Enrichment Analysis
The `enrichment` folder contains pipeline to calculate enrichment metrics like logAUC.

## Out-of-Sample Structure Predictions
The `out_of_sample` folder contains the workflow used to curate a dataset for structure prediction and then calculate relevant metrics.

## Prospective
The `prospective` folder contains example to run Boltz prediction and experimental data correlation calculation.

## Citing This Work
If you use this code or workflow in your research, please cite our paper:

> *Title*  
> *Authors*  
> *Journal / Year*
