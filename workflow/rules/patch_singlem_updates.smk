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
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time='08:00:00'
    shell:
        """
        #load conda environment
        module load conda/24.5.0
        source activate {config[conda_rairtable]}

        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/SMF/{config[version]} ; bye"
        rm -rf singlem_new.tsv

        for i in {config[workdir]}/output/*_readfraction.tsv;
            do sed '1d;' $i | sed 's/_M_1//' | cut -f1,4,5 >> singlem_new.tsv;
        done

        lftp sftp://erda -e "mirror -R {config[workdir]}/output /EarthHologenomeInitiative/Data/SMF/{config[version]}/; bye"

        Rscript {config[codedir]}/bonus/patch_singlem.R

        touch {output}
        """