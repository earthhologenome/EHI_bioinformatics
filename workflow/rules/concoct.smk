################################################################################
### Bin contigs using concoct
### Note that I've copied how the MetaWRAP pipeline calls MaxBin2 (https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-modules/binning.sh)
rule concoct:
    input:
        bam=os.path.join(
            config["workdir"], "bams/", "{PRB}_{EHI}_{EHA}.bam"
            ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{PRB}_{EHI}_{EHA}_binning/concoct_binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_binning"),
        contigsize=config["contigsize"]
    conda:
        f"{config['codedir']}/conda_envs/concoct.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/concoct_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/concoct_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with concoct"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt {params.contigsize} ]
        then
            touch {output}

        else
        
            for FILE in {input.bam}; do
                echo $FILE
                samtools index -@ {threads} -b $FILE
            done

            cut_up_fasta.py {input.contigs} -c 10000 --merge_last -b {params.outdir}/assembly_10K.bed -o 0 > {params.outdir}/assembly_10K.fa

            concoct_coverage_table.py {params.outdir}/assembly_10K.bed {input.bam} > {params.outdir}/concoct_depth.txt

            concoct \
            -l 1500 \
            -t {threads} \
            --coverage_file {params.outdir}/concoct_depth.txt \
            --composition_file {params.outdir}/assembly_10K.fa \
            -b {params.outdir}

        	merge_cutup_clustering.py {params.outdir}/clustering_gt1500.csv > {params.outdir}/clustering_gt1500_merged.csv
            mkdir -p {params.outdir}/concoct_bins 
            python {config[codedir]}/scripts/metawrap_split_concoct_bins.py {params.outdir}/clustering_gt1500_merged.csv {input.contigs} {params.outdir}/concoct_bins

            # Create output for the next rule
            touch {output}
        fi
        """