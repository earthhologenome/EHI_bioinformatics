################################################################################
### superPang
rule superpang:
    input:
        downloaded=os.path.join(
            config["magdir"],
            "mags_downloaded"
        )
    output:
        os.path.join(
            config["workdir"],
            "drep/",
            "figures/",
            config["dmb"] + "_Primary_clustering_dendrogram.pdf"
        )
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='08:00:00'
    benchmark:
        os.path.join(config["logdir"] + "/drep.benchmark.tsv")
    log:
        os.path.join(config["logdir"] + "/drep.log")
    message:
        "MESSAGE"
    shell:
        """
        # remove folder (in case of restart/server interruption)
        rm -rf {config[workdir]}/drep

        # Massage genome info file:
        sed 's/.fa/.fa.gz/g' mags.csv > mags_formatted.csv

        # Dereplicate these suckers:
        SuperPang.py \
                -f genome_path.tsv \
                -q genome_completeness.tsv \
                -t {threads} \
                -o output_dir \
                -u <header prefix>
                2> {log}

        for i in {config[workdir]}/drep/figures/*;
            do mv $i {config[workdir]}/drep/figures/{config[dmb]}_$(basename "$i");
        done

        #OUTPUTS: 
        - <name>_graph.fastg
        - <name>_graph_NBPorigins.csv

        #NEXT RULES = DRAM / mapping (NO GTDBTK) -> coverm

        #ALSO - R script for counts by contig && type [aux, core, sing] && genome <INPUT = COVERM OUTPUT>
        #ALSO API to pull GTDB taxonomy from AirTable

        tar -cvzf {config[workdir]}/drep/{config[dmb]}_drep_figures.tar.gz {config[workdir]}/drep/figures/

        """