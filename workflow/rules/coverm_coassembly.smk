################################################################################
### Calculate the number of reads that mapped to coassemblies
rule coverM_assembly:
    input:
        stats=os.path.join(
            config["workdir"],
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        contigmap=os.path.join(
            config["workdir"],
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.contigs",
            ),
        bam=os.path.join(
            config["workdir"], 
            "bams/", 
            "{PRB}_{EHI}_{EHA}.bam"
            ),
        contigs=os.path.join(
            config["workdir"], 
            "{EHA}_assembly/", 
            "{EHA}_contigs.fasta"
            )
    output:
        coverm=os.path.join(
            config["workdir"], 
            "coverm/", 
            "{PRB}_{EHI}_{EHA}_assembly_coverM.txt"
            ),
        euk=os.path.join(
            config["workdir"], 
            "coverm/", 
            "{PRB}_{EHI}_{EHA}_eukaryotic_coverM.tsv"
            ),
        tarball=os.path.join(
            config["workdir"], 
            "coverm/", 
            "{PRB}_{EHI}_{EHA}_coverm.tar.gz"
            )
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 4
    resources:
        load=1,
        mem_gb=64,
        time="00:15:00",
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


        # Tarball files then \:
        tar -cvzf {output.tarball} {output.coverm} {output.euk} {input.stats} {input.contigmap}
        """