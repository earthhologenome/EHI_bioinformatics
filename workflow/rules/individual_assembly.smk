################################################################################
### Perform individual assembly on each sample
rule assembly:
    input:
        r1=os.path.join(config["workdir"], "reads/", "{PRB}/", "{EHI}_M_1.fq.gz"),
        r2=os.path.join(config["workdir"], "reads/", "{PRB}/", "{EHI}_M_2.fq.gz"),
    output:
        os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta")
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=72,
        time=estimate_time_assembly
    benchmark:
        os.path.join(config["logdir"] + "/assembly_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/assembly_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Assembling {wildcards.EHA} using {params.assembler}"
    shell:
        """
        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -f \
                -o {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_assembly/
                2> {log}

        # Move the Coassembly to final destination
            mv {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_assembly/final.contigs.fa {output}

        # Reformat headers
            sed -i 's/ /-/g' {output}


        # if statement for situations where the assembly is too small to continue binning
        if [ $(( $(stat -c '%s' {output}) / 1024 / 1024 )) -lt 50 ]
        then
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_uploaded
            mkdir -p {config[workdir]}/bams
            touch {config[workdir]}/bams/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}.bam
            mkdir -p {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_binning
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_binning/binning_complete
            mkdir -p {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/
            touch {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv
            mkdir -p {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_refinement
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_refinement/{wildcards.EHA}_metawrap_50_10_bins.stats            

            # Create sample table for airtable
            echo -e "sample\tEHA_number\tEHI_number\tN50\tL50\tnum_contigs\tlargest_contig\tassembly_length\tnum_bins\tassembly_mapping_percent" > headers.tsv

            #Add sample IDs
            echo {wildcards.EHI}_{wildcards.EHA} >> {wildcards.EHA}_sample_ids.tsv
            echo {wildcards.EHA} >> {wildcards.EHA}_EHA_ids.tsv
            echo {wildcards.EHI} >> {wildcards.EHA}_EHI_ids.tsv

            echo -e "0\t0\t0\t0\t0\t0\t0" >> {wildcards.EHA}_null.tsv

            paste {wildcards.EHA}_sample_ids.tsv {wildcards.EHA}_EHA_ids.tsv {wildcards.EHA}_EHI_ids.tsv {wildcards.EHA}_null.tsv > {wildcards.EHA}_nullz.tsv
            cat headers.tsv {wildcards.EHA}_nullz.tsv > /projects/ehi/data/REP/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_final_stats.tsv


            ### Upload stats to AirTable:
            # python {config[codedir]}/airtable/add_asb_stats_airtable.py --report=/projects/ehi/data/REP/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_final_stats.tsv --asb={config[abb]}
        fi

        """