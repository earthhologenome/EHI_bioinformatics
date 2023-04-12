################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        bam=os.path.join(config["workdir"], "{PRB}", "{EHI}", "{EHI}", "{EHA}.bam"),
        contigs=os.path.join(config["workdir"], "{PRB}" "{EHI}" "{EHA}_contigs.fasta"),
    output:
        os.path.join(
            config["workdir"], "{PRB}", "{EHI}", "{EHA}_binning/binning_complete"
        ),
    params:
        outdir=directory("{config['workdir']}/{PRB}/{EHI}/{EHA}_binning"),
    conda:
        f"{config['codedir']}/conda_envs/metawrap.yaml"
    threads: 16
    resources:
        mem_gb=96,
        time="06:00:00",
    benchmark:
        "{{config['logdir']}}/binning_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/binning_log_{PRB}_{EHI}_{EHA}.log",
    message:
        "Binning {wildcards.EHA} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
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