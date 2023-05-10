################################################################################
### Create summary table from outputs
rule report:
    input:
        coverm=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}_coverM_mapped_host.tsv"
            ),
            sample=SAMPLE
        ),
        fastp=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}.json"
            ),
            sample=SAMPLE
        ),
        read_fraction=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}_readfraction.tsv"
            ),
            sample=SAMPLE
        ),        
        uploaded=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}_uploaded"
            ),
            sample=SAMPLE
        ),
        npstats=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}_np.tsv"
            ),
            sample=SAMPLE
        )
    output:
        report=os.path.join(
            "/projects/ehi/data/REP/",
            config["prb"] + ".tsv"
        ),
        npar_metadata=os.path.join(
            config["workdir"],
            config["prb"] + "_nonpareil_metadata.tsv"
        )
    params:
        tmpdir=os.path.join(
            config["workdir"],
            "tmp/"
        ),
        npar=expand(
            os.path.join(
                config["workdir"],
                "misc/{sample}.npo"
            ),
            sample=SAMPLE
        ),
        misc_dir=os.path.join(
            config["workdir"],
            "misc/"
        )
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=1,
        mem_gb=24,
        time='00:20:00'
    message:
        "Creating a final preprocessing report"
    shell:
        """
        #Create nonpareil sample metadata file
        mkdir -p {params.tmpdir}
        for i in {params.npar}; do echo $(basename $i) >> {params.tmpdir}/files.txt; done
        for i in {params.npar}; do echo $(basename ${{i/.npo/}}) >> {params.tmpdir}/names.txt; done
        for i in {params.npar}; do echo "#f03b20" >> {params.tmpdir}/colours.txt; done
        echo -e "File\tName\tColour" > {params.tmpdir}/headers.txt
        paste {params.tmpdir}/files.txt {params.tmpdir}/names.txt {params.tmpdir}/colours.txt > {params.tmpdir}/merged.tsv
        cat {params.tmpdir}/headers.txt {params.tmpdir}/merged.tsv > {output.npar_metadata}

        #Create preprocessing report
        mkdir -p {params.tmpdir}
        for i in {input.coverm}; do echo $(basename ${{i/_coverM_mapped_host.tsv}}) >> {params.tmpdir}/names.tsv; done
        for i in {input.coverm}; do grep -v 'Genome' $i | grep -v 'unmapped' | cut -f3; done >> {params.tmpdir}/host_reads.tsv

        for i in {input.fastp}; do grep '"total_reads"' $i | sed -n 1p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/read_pre_filt.tsv
        for i in {input.fastp}; do grep '"total_reads"' $i | sed -n 2p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/read_post_filt.tsv
        for i in {input.fastp}; do grep '"total_bases"' $i | sed -n 1p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/bases_pre_filt.tsv
        for i in {input.fastp}; do grep '"total_bases"' $i | sed -n 2p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/bases_post_filt.tsv
        for i in {input.fastp}; do grep 'adapter_trimmed_reads' $i | cut -f2 --delimiter=: | tr -d ',' | tr -d ' '; done >> {params.tmpdir}/adapter_trimmed_reads.tsv
        for i in {input.fastp}; do grep 'adapter_trimmed_bases' $i | cut -f2 --delimiter=: | tr -d ',' | tr -d ' '; done >> {params.tmpdir}/adapter_trimmed_bases.tsv

        #parse singlem estimates
        for i in {input.read_fraction}; do sed '1d;' $i | cut -f2,3,4 >> {params.tmpdir}/singlem.tsv; done

        #parse nonpareil estimates
        for i in {input.npstats}; do sed '1d;' $i | cut -f2,3,4,5,6,7 >> {params.tmpdir}/npstats.tsv; done

        paste {params.tmpdir}/names.tsv {params.tmpdir}/read_pre_filt.tsv {params.tmpdir}/read_post_filt.tsv {params.tmpdir}/bases_pre_filt.tsv {params.tmpdir}/bases_post_filt.tsv {params.tmpdir}/adapter_trimmed_reads.tsv {params.tmpdir}/adapter_trimmed_bases.tsv {params.tmpdir}/host_reads.tsv {params.tmpdir}/singlem.tsv {params.tmpdir}/npstats.tsv > {params.tmpdir}/preprocessing_stats.tsv
        echo -e "EHI_number\treads_pre_fastp\treads_post_fastp\tbases_pre_fastp\tbases_post_fastp\tadapter_trimmed_reads\tadapter_trimmed_bases\thost_reads\tbacterial_archaeal_bases\tmetagenomic_bases\tsinglem_fraction\tkappa\tC\tLR\tmodelR\tLRstar\tdiversity" > {params.tmpdir}/headers.tsv
        cat {params.tmpdir}/headers.tsv {params.tmpdir}/preprocessing_stats.tsv > {output.report}

        cp {output.report} {params.misc_dir}
        cp {output.npar_metadata} {params.misc_dir}
        tar -czf {config[workdir]}/{config[prb]}_stats.tar.gz {params.misc_dir}

        #Upload stats and report to ERDA for storage
        lftp sftp://erda -e "put {config[workdir]}/{config[prb]}_stats.tar.gz -o /EarthHologenomeInitiative/Data/PPR/{config[prb]}/; bye"
        sleep 10
        lftp sftp://erda -e "put {output.report} -o /EarthHologenomeInitiative/Data/REP/; bye"

        #Automatically update the AirTable with the preprocessing stats
        python {config[codedir]}/airtable/add_prb_stats_airtable.py --report={output.report} --prb={config[prb]} 

        #Indicate that the PRB is done in AirTable
        python {config[codedir]}/airtable/log_prb_done_airtable.py --code={config[prb]}
       
        #Clean up the files/directories
        rm {config[workdir]}/{config[prb]}_stats.tar.gz
        rm -r {config[workdir]}/{config[hostgenome]}/
        rm -r {params.misc_dir}
        rm -r {params.tmpdir}

        """