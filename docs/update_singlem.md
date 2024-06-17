# Guide for updating SingleM values in AirTable
Raphael Eisenhofer 2024/06

### First, create a new conda environment with the latest version of SingleM
```
conda create --prefix /projects/ehi/data/SMF/conda/singlem_0.18.0 singlem=0.18.0
```

### Enter the environment and download the latest metapackage
```
conda activate /projects/ehi/data/SMF/conda/singlem_0.18.0

singlem data --output-directory /projects/ehi/data/SMF/metapackages/

```

### Load the rairtable conda environment and get inputs
```
conda activate /projects/ehi/data/SMF/conda/rairtable

Rscript get_singlem_stats_and_input.R

conda deactivate
```

### Execute the snakefile
```

```