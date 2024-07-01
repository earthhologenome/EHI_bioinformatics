################################################################################
### Register samples and upload reads to ENA (also push ENA codes to airtable)
rule register_sample_ena_upload_reads:
    input:
        sample_checklist=os.path.join(
            config["workdir"],
            "{EHI}_sample_checklist.tsv"
        ),
        experiment_checklist=os.path.join(
            config["workdir"],
            "{EHI}_experiment_checklist.tsv"
        ),
        run_checklist=os.path.join(
            config["workdir"],
            "{EHI}_run_checklist.tsv"
        )
    output:
        sample_checklist_updated=os.path.join(
            config["workdir"],
            "{EHI}_sample_checklist_updated.tsv"
        ),
        experiment_checklist_updated=os.path.join(
            config["workdir"],
            "{EHI}_experiment_checklist_updated.tsv"
        ),
        run_checklist_updated=os.path.join(
            config["workdir"],
            "{EHI}_run_checklist_updated.tsv"
        ),
        accessions_uploaded=os.path.join(
            config["workdir"],
            "{EHI}_accessions_uploaded"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time=estimate_time_download
    message:
        "Regsitering {wildcards.EHI} to ENA and uploading raw reads"
    shell:
        """
        source activate /projects/ehi/data/0_Environments/conda/ena_upload

        #Register the samples and upload the reads to the ENA
        ena-upload-cli \
        --action add \
        --center 'Earth Hologenome Initiative' \
        --sample {input.sample_checklist} \
        --experiment {input.experiment_checklist} \
        --run {input.run_checklist} \
        --checklist ERC000013 \
        --data {config[workdir]}/{wildcards.EHI}*.fq.gz \
        --secret /projects/ehi/data/.secret.yml

        #Use API to patch the ENA sample accession to the EHI AirTable (Samples table)
        python {config[codedir]}/airtable/add_ena_sample_accession.py \
        --sample `sed '1d;' {output.sample_checklist_updated} | cut -f1` \
        --sample_acc `sed '1d;' {output.sample_checklist_updated} | cut -f20`

        #Use API to patch the ENA experiment and run accessions to the EHI AirTable ('SE Samples' table)
        python {config[codedir]}/airtable/add_ena_sample_accession.py \
        --ehi `sed '1d;' {output.experiment_checklist_updated} | cut -f2` \
        --exp_acc `sed '1d;' {output.experiment_checklist_updated} | cut -f16` \
        --run_acc `tail -1 {output.run_checklist_updated} | cut -f5 `

        #Close job
        touch {output.accessions_uploaded}
        """