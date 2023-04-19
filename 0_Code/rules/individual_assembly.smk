################################################################################
### Perform individual assembly on each sample
rule assembly:
    input:
        r1=os.path.join(config["workdir"], "{PRB}/", "{EHI}_M_1.fq.gz"),
        r2=os.path.join(config["workdir"], "{PRB}/", "{EHI}_M_2.fq.gz"),
    output:
        os.path.join(config["workdir"], "{PRB}_{EHA}_assembly", "{EHA}_contigs.fasta")
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=128,
        time="12:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/assembly_benchmark_{PRB}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/assembly_log_{PRB}_{EHA}.log")
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
                -o {config[workdir]}/{wildcards.EHA}/assembly/
                2> {log}

        # Move the Coassembly to final destination
            mv {config[workdir]}/{wildcards.EHA}/assembly/final.contigs.fa {output}

        # Reformat headers
            sed -i 's/ /-/g' {output}
        """