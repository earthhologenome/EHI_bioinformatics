################################################################################
### Fetch MAGs from ERDA and mag checklist from AirTable
rule download_mags:
    input:
        accessions_uploaded=os.path.join(
            config["workdir"],
            "{EHI}_accessions_uploaded"
        ),
        sample_checklist_updated=os.path.join(
            config["workdir"],
            "{EHI}_sample_checklist_updated.tsv"
        ),
    output:
        mag_checklist=os.path.join(
            config["workdir"],
            "{EHI}_mag_checklist.tsv"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time=estimate_time_download
    message:
        "Fetching {wildcards.EHI} MAGs ERDA and creating MAG checklist"
    shell:
        """
        #Fetch most of the checklist from AirTable
        python {config[codedir]}/airtable/get_mags_ena.py \
        --ehi {wildcards.EHI} \
        --ehs {config[ehs]}

        #Get the rest of the checklist metadata from the previously generated sample checklist
        {input.sample_checklist_updated}


        """