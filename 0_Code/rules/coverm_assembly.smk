################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        stats=os.path.join(
            config["workdir"],
            "{PRB}",
            "{EHI}",
            "{EHA}_refinement",
            "{EHA}_metawrap_70_10_bins.stats",
        ),
        contigmap=os.path.join(
            config["workdir"],
            "{PRB}",
            "{EHI}",
            "{EHA}_refinement",
            "{EHA}_metawrap_70_10_bins.contigs",
        ),
        bam=os.path.join(config["workdir"], "{PRB}", "{EHI}", "{EHI}", "{EHA}.bam"),
        contigs=os.path.join(config["workdir"], "{PRB}" "{EHI}" "{EHA}_contigs.fasta"),
    output:
        coverm=os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_assembly_coverM.txt"
        ),
        euk=os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_eukaryotic_coverM.tsv"
        ),
        contigs_gz=os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_contigs.fasta.gz"
        ),
    params:
        refinement_files="{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/",
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 8
    resources:
        mem_gb=64,
        time="00:30:00",
    benchmark:
        "{{config['logdir']}}/coverm_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/coverm_log_{PRB}_{EHI}_{EHA}.log",
    message:
        "Calculating assembly mapping rate for {wildcards.EHA} with CoverM"
    shell:
        """
        coverm genome \
            -b {input.bam} \
            --genome-fasta-files {input.contigs} \
            -m relative_abundance \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.coverm}

        #Run coverm for the eukaryotic assessment pipeline
        coverm genome \
            -s - \
            -b {input.bam} \
            -m relative_abundance count mean covered_fraction \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output.euk}

        # Compress the contigs
        pigz -p {threads} {input.contigs}

        #Print the number of MAGs to a file for combining with the assembly report
        ls -l {params.refinement_files}/metawrap_70_10_bins/*.fa.gz | wc -l > {wildcards.EHA}_bins.tsv;
        """