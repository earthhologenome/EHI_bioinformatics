################################################################################
### Perform individual assembly on each sample
rule assembly:
    input:
        r1=os.path.join(config["workdir"], "{PRB}/", "{EHI}_M_1.fq.gz"),
        r2=os.path.join(config["workdir"], "{PRB}/", "{EHI}_M_2.fq.gz"),
    output:
        os.path.join(config["workdir"], "{PRB}/" "{EHI}/" "{EHA}_contigs.fasta"),
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=128,
        time="01:00:00",
    benchmark:
        "{{config['logdir']}}/assembly_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/assembly_log_{PRB}_{EHI}_{EHA}.log",
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
                -o {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}
                2> {log}

        # Move the Coassembly to final destination
            mv {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/final.contigs.fa {output}

        # Reformat headers
            sed -i 's/ /-/g' {output}
        """