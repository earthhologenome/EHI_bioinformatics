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
snakemake \
    -s /projects/ehi/data/0_Code/EHI_bioinformatics_1.1/bonus/update_singlem_airtable.snakefile \
    --configfile /projects/ehi/data/0_Code/EHI_bioinformatics_1.1/config/singlem_update.yaml \
    -j 2 \
    --cluster "sbatch --mem {resources.mem_gb}G -c {threads} --time {resources.time} -v" \
    --use-conda \
    --conda-frontend conda \
    --latency-wait 600 
```