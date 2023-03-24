# EHI_Bioinformatics 
# 🐨->💩->🦠->🧬->🖥️->😏
Bioinformatics pipeline to process EHI data.

*updated 24/03/2023, Raphael Eisenhofer*

#### General information:
This pipeline uses [![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io), and manages dependencies using conda (or mamba) for reproducibility and deployability. The 0_Code directory contains the snakefiles, scripts, and conda environment yamls. 

#### Getting started:
Firstly, you'll need to set up an alias for connecting to ERDA -- this is **essential**.