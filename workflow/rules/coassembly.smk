################################################################################
### Perform coassembly on samples
rule coassembly:
    input:
        r1=expand(
            os.path.join(
                config["workdir"],
                "reads/",
                "{combo[0]}/",
                "{combo[1]}_M_1.fq.gz"
                ),
                combo=valid_combinations
        ),
        r2=expand(
            os.path.join(
                config["workdir"],
                "reads/",
                "{combo[0]}/",
                "{combo[1]}_M_2.fq.gz"
                ),
                combo=valid_combinations
        )   
    output:
        os.path.join(config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta")
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=128,
        time="16:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/assembly_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/assembly_log_{EHA}.log")
    message:
        "Assembling {wildcards.EHA} using {params.assembler}"
    shell:
        """
        # Setup input variables for megahit
        R1=$(for i in {input.r1}; do echo $i | tr '\n' ,; done)
        R2=$(for i in {input.r2}; do echo $i | tr '\n' ,; done)

        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 $R1 -2 $R2 \
                -f \
                -o {config[workdir]}/{wildcards.EHA}_assembly/
                2> {log}

        # Move the Coassembly to final destination
        mv {config[workdir]}/{wildcards.EHA}_assembly/final.contigs.fa {output}

        # Reformat headers
        sed -i 's/ /-/g' {output}
        """