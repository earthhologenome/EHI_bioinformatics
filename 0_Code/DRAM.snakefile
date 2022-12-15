rule all:
    input:
        "3_Outputs/test_bin/annotations.tsv")

################################################################################
### Functionally annotate MAGs with DRAM
rule DRAM:
    input:
        bin = "3_Outputs/test_bin.fa.gz"
    output:
        annotation = "3_Outputs/test_bin/annotations.tsv"),
    params:
        database = expand("{database}", database=config['database']),
    conda:
        "conda_envs/DRAM.yaml"
    threads:
        24
    resources:
        mem_gb=64,
        time='03:00:00'
    benchmark:
        "3_Outputs/0_Logs/DRAM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/DRAM.log"
    message:
        "Using DRAM to functionally annotate bin"
    shell:
        """
        DRAM-setup.py import_config --config_loc /projects/mjolnir1/people/ncl550/0_software/20210705.dram.config

        DRAM.py annotate \
            -i {input.bin} \
            -o test_bin \
            --threads {threads} \
            --min_contig_size 1500 

        """