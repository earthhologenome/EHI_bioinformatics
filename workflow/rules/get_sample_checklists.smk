################################################################################
### Get ENA checklist information for sample
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
        
        sed -i'' "s/'//g" {output.sample_checklist}
        sed -i'' "s/\[//g" {output.sample_checklist}
        sed -i'' "s/\]//g" {output.sample_checklist}

        #pull experiment checklist
        python {config[codedir]}/airtable/get_experiment_checklist_airtable.py \
        --ehi={wildcards.EHI}

        sed -i'' "s/'//g" {output.experiment_checklist}
        sed -i'' "s/\[//g" {output.experiment_checklist}
        sed -i'' "s/\]//g" {output.experiment_checklist}

        #make run checklist
        echo -e "alias\texperiment_alias\tfile_name\tfile_type" >> run_header.tsv
        echo -e "{wildcards.ehi}\t`grep {wildcards.EHI} ehi_numbers.tsv | cut -f2`\t{input.r1}\tfastq" > run_r1.tsv
        echo -e "{wildcards.ehi}\t`grep {wildcards.EHI} ehi_numbers.tsv | cut -f2`\t{input.r2}\tfastq" > run_r2.tsv
        cat run_header.tsv run_r1.tsv run_r2.tsv > {output.run_checklist}

        """