################################################################################
### Combine contigs from different individual assemblies
rule combine_contigs:
    input:
        expand(os.path.join(
            config["workdir"], 
            "{combo[0]}_{combo[1]}_assembly/",
            "{combo[2]}_contigs.fasta"
                ), combo=valid_combinations
            ),
    output:
        os.path.join(config["workdir"], "{EHA}_combined_contigs.fasta.gz")
    params:
        assembler=expand("{assembler}", assembler=config["assembler"]),
    # conda:
    #     f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 1
    resources:
        mem_gb=16,
        time='00:30:00'
    message:
        "Combining contigs into a single fasta file"
    shell:
        """
        source activate /projects/mjolnir1/people/ncl550/0_software/miniconda/envs/semibin_1.5.1

        SemiBin2 concatenate_fasta \
        -i {input} \
        -o {output} \
        -s : \
        -m 2000
        
        """