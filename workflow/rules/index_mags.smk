################################################################################
### Index mags
rule mag_index:
    input:
        expand(
            os.path.join(
                config["magdir"], 
                "{mag}"
            ), mag = MAG
        )
    output:
        mags=os.path.join(
            config["magdir"],
            config["dmb"] + "_mags.fasta.gz"
        ),
        index=os.path.join(
            config["magdir"],
            config["dmb"] + "_mags.fasta.gz.rev.2.bt2l"
        )
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 16
    resources:
        mem_gb=96,
        time="02:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/index_mags_benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/index_mags_log.log")
    message:
        "Indexing mags using Bowtie2"
    shell:
        """
        # Reformat MAG headers for CoverM
        for mag in {input};
            do rename.sh \
                in=$mag \
                out=${{mag/.fa.gz/_renamed.fa.gz}} \
                zl=9 \
                prefix=$(basename ${{mag/.fa.gz/^}});
        done

        # Concatenate MAGs
        cat {config[magdir]}/*_renamed.fa.gz > {config[magdir]}/{config[dmb]}_mags.fasta.gz

        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.mags} {output.mags} \
        &> {log}
        """