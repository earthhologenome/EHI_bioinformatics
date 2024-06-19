################################################################################
### Fetch preprocessed reads from ERDA
rule download_reads_singlem:
    output:
        r1=os.path.join(
            config["workdir"], 
            "reads/", 
            "{EHI}_M_1.fq.gz"
        ),
        r2=os.path.join(
            config["workdir"], 
            "reads/", 
            "{EHI}_M_2.fq.gz"
        )
    threads: 1
    resources:
        load=8,
        mem_gb=8,
        time='01:00:00'
    message:
        "Fetching metagenomics reads for {wildcards.EHI} from ERDA"
    shell:
        """
        wget --no-verbose `grep '{wildcards.EHI}' singlem_input.csv | cut -f4 -d , | sed 's/"//g'`
        mv {wildcards.EHI}_M_1.fq.gz {output.r1}

        wget --no-verbose `grep '{wildcards.EHI}' singlem_input.csv | cut -f5 -d , | sed 's/"//g'`
        mv {wildcards.EHI}_M_2.fq.gz {output.r2}

        """