################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        bam=expand(os.path.join(
            config["workdir"], "bams/", "{combo[0]}_{combo[1]}_{combo[2]}.bam"
            ),
            combo=valid_combinations
        ),
        contigs=os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{EHA}_binning/binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{EHA}_binning")
    threads: 16
    resources:
        mem_gb=96,
        time="16:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/binning_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/binning_log_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        #Installing metawrap via conda is a pain in the arse, so using the module on Mjolnir here.
        #It should be possible to get the conda environment done manually 'mamba create -n X -c ursky metawrap-mg'
        module load metawrap-mg/1.3.2
        #seems to be an issue with mamba=1.4.1 module, so unload it to get concoct to work ^_^"
        module unload mamba

        # Create dummy fq/assembly files to trick metaWRAP into running without mapping
        mkdir -p {params.outdir}
        mkdir -p {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input.bam}; do ln -sf $bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {resources.mem_gb} \
            -a {input.contigs} \
            -l 1500 \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq

        # Create output for the next rule
        touch {output}
        """