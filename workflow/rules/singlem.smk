################################################################################
### Estimate the fraction of bacterial and archaeal DNA using SingleM read fraction
rule singlem:
    input:
        npo=os.path.join(
            config["workdir"],
            "misc/{sample}.npo"
        ),
        non_host_r1=os.path.join(
            config["workdir"],
            "{sample}_M_1.fq.gz"
        ),
        non_host_r2=os.path.join(
            config["workdir"],
            "{sample}_M_2.fq.gz"
        )
    output:
        otu=os.path.join(
            config["workdir"],
            "misc/{sample}_OTU.tsv.gz"
        ),
        condense=os.path.join(
            config["workdir"],
            "misc/{sample}_condense.tsv"
        ),
        read_fraction=os.path.join(
            config["workdir"],
            "misc/{sample}_readfraction.tsv"
        )
    params:
        otu_uncompressed=os.path.join(
            config["workdir"],
            "misc/{sample}_OTU.tsv"
        ),
        read_fraction_taxa=os.path.join(
            config["workdir"],
            "misc/{sample}_readfraction_per_taxa.tsv"
        ),
        archive=os.path.join(
            config["workdir"],
            "misc/{sample}_archive.json"
        ),
# Current issue with snakemake and pre-built conda environments: https://github.com/snakemake/snakemake/pull/1708
    conda:
        f"{config['codedir']}/conda_envs/singlem.yaml"
    threads:
        3
    resources:
        load=1,
        mem_gb=8,
        time=estimate_time_singlem
    benchmark:
        os.path.join(config["logdir"] + "/{sample}_singlem.benchmark.tsv")
    message:
        "Estimating microbial fraction using singlem"
    shell:
        """
        #Temp fix until snakemake is fixed or singlem conda recipe is updated
        export PATH='/projects/ehi/data/0_Environments/github_repos/singlem/bin':$PATH
        export SINGLEM_METAPACKAGE_PATH='/projects/ehi/data/0_Environments/databases/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'

        #Try to fix /tmp folder running out of space:
        export TMPDIR={config[workdir]}/tmpdir
        mkdir -p $TMPDIR

        #IF statement to account for situations where there are not enough
        #microbial reads in a sample (e.g. high host% or non-metagenomic sample)
        #In this case, if R1 has > 150 Mbytes, run, else, skip:

        if [ $(( $(stat -c '%s' {input.non_host_r1}) / 1024 / 1024 )) -gt 150 ]
        then

        #Run singlem pipe
        singlem pipe \
            -1 {input.non_host_r1} \
            -2 {input.non_host_r2} \
            --otu-table {params.otu_uncompressed} \
            --taxonomic-profile {output.condense} \
            --archive-otu-table {params.archive} \
            --threads {threads}

        #Compress pipe file
        gzip -f {params.otu_uncompressed}
        gzip -f {params.archive}

            #IF statement for files without data
            if [ $(( $(stat -c '%s' {output.condense}) )) -eq 25 ]
            then
            echo -e "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\taverage_genome_size\n0\t0\t0\t0.0\t0" > {output.read_fraction}
            
            else        
            #Run singlem read_fraction
            singlem microbial_fraction \
                -1 {input.non_host_r1} \
                -2 {input.non_host_r2} \
                --input-profile {output.condense} \
                --output-tsv {output.read_fraction} \
                --output-per-taxon-read-fractions {params.read_fraction_taxa}

            #And upload untarred files for easy access
            sftp erda:/EarthHologenomeInitiative/Data/PPR/{config[prb]} <<< $'put {config[workdir]}/misc/{wildcards.sample}*'

            fi

        #Otheriwse, don't run singlem
        else
        echo "SingleM analysis not performed"
        touch {output.condense}
        touch {output.otu}
        echo -e "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\taverage_genome_size\n0\t0\t0\t0.0\t0" > {output.read_fraction}
        
        fi

        #If statement for cases when singlem does not produce a condense output
        if [ -f {params.read_fraction_taxa} ]
        then
        #Compress read_fraction_per_taxa file
        gzip -f {params.read_fraction_taxa}

        else
        echo "no microbes in sample"
        fi
        
        """