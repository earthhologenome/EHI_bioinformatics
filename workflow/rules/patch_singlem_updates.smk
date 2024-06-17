###############################################################################
## Upload SingleM outputs and patch new values to AirTable
rule patch_singlem_updates:
    input:
        singlem=expand(
            os.path.join(
                config["workdir"], "output/{EHI}_readfraction.tsv"
                ),
                EHI=EHI
        )
    output:
        os.path.join(
            config["workdir"],
            "singlem_updated"
        )
    params:
        stats_dir=directory(os.path.join(
            config["workdir"],
            "stats/")
            ),
    conda:
        f"{config['conda_rairtable']}"
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time='08:00:00'
    shell:
        """
        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/SMF/{config[version]} ; bye"
        rm -rf {params.stats_dir}
        mkdir -p {params.stats_dir}

        for i in {config[workdir]}/output/*_readfraction.tsv;
            do XXX >> {params.stats_dir}/annos.tsv;
        done

        lftp sftp://erda -e "put {params.stats_dir} -o /EarthHologenomeInitiative/Data/SMF/{config[version]}/; bye"

        lftp sftp://erda -e "put outputs/ -o /EarthHologenomeInitiative/Data/SMF/{config[version]}/; bye"

        Rscript {config[codedir]}/bonus/patch_singlem.R

        touch {output}
        """