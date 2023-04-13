################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        stats=os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
        ),
        contigmap=os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.contigs",
        ),
        bam=os.path.join(config["workdir"], "{PRB}/", "{EHI}/", "{EHI}_{EHA}.bam"),
        contigs=os.path.join(config["workdir"], "{PRB}/" "{EHI}/" "{EHA}_contigs.fasta"),
    output:
        coverm=os.path.join(
            config["workdir"], "{PRB}/", "{EHI}/", "{EHA}_assembly_coverM.txt"
        ),
        euk=os.path.join(
            config["workdir"], "{PRB}/", "{EHI}/", "{EHA}_eukaryotic_coverM.tsv"
        ),
        contigs_gz=os.path.join(
            config["workdir"], "{PRB}/", "{EHI}/", "{EHA}_contigs.fasta.gz"
        ),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 8
    resources:
        load=8,
        mem_gb=64,
        time="00:30:00",
    benchmark:
        os.path.join(config["logdir"] + "/coverm_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/coverm_log_{PRB}_{EHI}_{EHA}.log")
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

        # Tarball files then upload to ERDA:
        tar -cvzf {wildcards.EHA}_coverm.tar.gz {output.coverm} {output.euk}
        lftp sftp://erda -e "put {wildcards.EHA}_coverm.tar.gz -o /EarthHologenomeInitiative/Data/ASB/{config[abb]}/; bye"
        """