# Guide for updating SingleM values in AirTable
Raphael Eisenhofer 2024/06

### Intro
Use this guide if you want to update SingleM and the associated metapackages.

SingleM outputs and stats are stored on ERDA:
`/EarthHologenomeInitiative/Data/SMF/`


### First, create a new conda environment with the latest version of SingleM
```
mamba create --prefix /projects/ehi/data/SMF/conda/singlem_0.18.0 singlem=0.18.0
```

### Enter the environment and download the latest metapackage
```
conda activate /projects/ehi/data/SMF/conda/singlem_0.18.0

singlem data --output-directory /projects/ehi/data/SMF/metapackages/

#set path to metapackage variable in conda
conda env config vars set SINGLEM_METAPACKAGE_PATH=/projects/ehi/data/SMF/metapackages/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb

```

### Load the rairtable conda environment and get inputs
```
cd /projects/ehi/data/SMF

conda activate /projects/ehi/data/SMF/conda/rairtable

Rscript /projects/ehi/data/0_Code/EHI_bioinformatics_1.1/bonus/get_singlem_stats_and_input.R

conda deactivate
```

### Execute the snakefile
Had to edit conda.py in /projects/ehi/data/SMF/conda/snakemake, see https://github.com/snakemake/snakemake/pull/1708/files
```
conda activate /projects/ehi/data/SMF/conda/snakemake

cd /projects/ehi/data/SMF

snakemake \
    -s /projects/ehi/data/0_Code/EHI_bioinformatics_1.1/bonus/update_singlem_airtable.snakefile \
    --configfile /projects/ehi/data/0_Code/EHI_bioinformatics_1.1/config/singlem_update.yaml \
    -j 100 \
    --cluster "sbatch --mem {resources.mem_gb}G -c {threads} --time {resources.time} --job-name=EHI-SMF-{rule} --parsable -v" \
    --retries 3 \
    --latency-wait 600 
```

# Patch AirTable
Note that for some reason I couldn't get this to work on mjolnir, I was getting this error:
`Error: Error in PATCH Error Code 422: Unprocessable entity. The request was well-formed but was unable to be followed due to semantic errors.`

Therefore, you'll need to download the output data locally and run the following Rscript:

```


```