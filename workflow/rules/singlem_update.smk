################################################################################
### Estimate the fraction of bacterial and archaeal DNA using SingleM read fraction
rule singlem:
    input:
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
    output:
        archive_gz=os.path.join(
            config["workdir"],
            "output/{EHI}_archive.json.gz"
        ),
        profile=os.path.join(
            config["workdir"],
            "output/{EHI}_profile.tsv"
        ),
        otu_table=os.path.join(
            config["workdir"],
            "output/{EHI}_OTU_table.tsv"
        ),
        read_fraction=os.path.join(
            config["workdir"],
            "output/{EHI}_readfraction.tsv"
        ),
        per_taxon=os.path.join(
            config["workdir"],
            "output/{EHI}_pertaxon.tsv"
        )
    params:
        archive=os.path.join(
            config["workdir"],
            "output/{EHI}_archive.json"
        ),
    conda:
        f"{config['conda_singlem']}"
    threads:
        2
    resources:
        load=1,
        mem_gb=8,
        time=estimate_time_singlem
    message:
        "Estimating microbial fraction using singlem"
    shell:
        """
        #Try to fix /tmp folder running out of space:
        export TMPDIR={config[workdir]}/tmpdir
        mkdir -p $TMPDIR

        #Run singlem pipe
        singlem pipe \
            -1 {input.r1} \
            -2 {input.r2} \
            --otu-table {output.otu_table} \
            --taxonomic-profile {output.profile} \
            --archive-otu-table {params.archive} \
            --threads {threads}

        #Run singlem microbial_fraction
        singlem microbial_fraction \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {output.profile} \
            --output-tsv {output.read_fraction} \
            --output-per-taxon-read-fractions {output.per_taxon}
        

        #compress some files
        gzip {params.archive}

        """