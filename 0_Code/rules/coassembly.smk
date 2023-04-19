################################################################################
### Perform coassembly on samples
rule coassembly:
    input:
        r1=expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}_M_1.fq.gz"
                ),
                combo=valid_combinations
        ),
        r2=expand(
            os.path.join(
                config["workdir"],
                "{combo[0]}/",
                "{combo[1]}_M_2.fq.gz"
                ),
                combo=valid_combinations
        )   
    output:
        os.path.join(config["workdir"], "{EHA}_contigs.fasta")
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=128,
        time="12:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/assembly_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/assembly_log_{EHA}.log")
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
        mv {config[workdir]}/final.contigs.fa {output}

        # Reformat headers
        sed -i 's/ /-/g' {output}
        """