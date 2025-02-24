################################################################################
### Get ENA checklist information for sample
## Note, I added an '_ena' prefix to the EHI numbers in the experiment/run checklists.
## This is to fix an issue with a previous variant of the ENA pipeline already creating accessions -
## with these codes. 
rule get_sample_checklists:
    input:
        r1=os.path.join(
            config["workdir"],
            "{EHI}_raw_1.fq.gz"
        ),
        r2=os.path.join(
            config["workdir"],
            "{EHI}_raw_2.fq.gz"
        )
    output:
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
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time="00:10:00"
    message:
        "Fetching sample checklist for {wildcards.EHI}"
    shell:
        """
        #pull sample checklist
        python {config[codedir]}/airtable/get_sample_checklist_airtable.py \
        --ehi={wildcards.EHI} \
        --sample=`grep {wildcards.EHI} ehi_numbers.tsv | cut -f2`

        mv {wildcards.EHI}_sample_checklist.tsv {output.sample_checklist}
        
        sed -i'' "s/'//g" {output.sample_checklist}
        sed -i'' "s/\[//g" {output.sample_checklist}
        sed -i'' "s/\]//g" {output.sample_checklist}
        sed -i'' 's/TEMP/{wildcards.EHI}/g' {output.sample_checklist}

        #pull experiment checklist
        python {config[codedir]}/airtable/get_experiment_checklist_airtable.py \
        --ehi={wildcards.EHI}

        mv {wildcards.EHI}_experiment_checklist.tsv {output.experiment_checklist}

        sed -i'' "s/'//g" {output.experiment_checklist}
        sed -i'' "s/\[//g" {output.experiment_checklist}
        sed -i'' "s/\]//g" {output.experiment_checklist}
        sed -i'' 's/EHI/ena_EHI/2;s/EHI/ena_EHI/' {output.experiment_checklist}

        #make run checklist
        echo -e "alias\texperiment_alias\tfile_name\tfile_type" > {wildcards.EHI}_run_header.tsv
        echo -e "ena_{wildcards.EHI}\t`echo ena_{wildcards.EHI}`\t{wildcards.EHI}_raw_1.fq.gz\tfastq" > {wildcards.EHI}_run_r1.tsv
        echo -e "ena_{wildcards.EHI}\t`echo ena_{wildcards.EHI}`\t{wildcards.EHI}_raw_2.fq.gz\tfastq" > {wildcards.EHI}_run_r2.tsv
        cat {wildcards.EHI}_run_header.tsv {wildcards.EHI}_run_r1.tsv {wildcards.EHI}_run_r2.tsv > {output.run_checklist}

        """