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
        pipe=os.path.join(
            config["workdir"],
            "misc/{sample}_pipe.tsv.gz"
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
        pipe_uncompressed=os.path.join(
            config["workdir"],
            "misc/{sample}_pipe.tsv"
        ),
        read_fraction_taxa=os.path.join(
            config["workdir"],
            "misc/{sample}_readfraction_per_taxa.tsv"
        )
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
        export SINGLEM_METAPACKAGE_PATH='/projects/ehi/data/0_Environments/databases/S3.1.0.metapackage_20221209.smpkg.zb/'

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
            --otu-table {params.pipe_uncompressed} \
            --taxonomic-profile {output.condense} \
            --threads {threads}

        #Compress pipe file
        gzip {params.pipe_uncompressed}

            #IF statement for files without data
            if [ $(( $(stat -c '%s' {output.condense}) )) -eq 25 ]
            then
            echo -e "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\n0\t0\t0\t0.0%" > {output.read_fraction}
            
            else        
            #Run singlem read_fraction
            singlem read_fraction \
                -1 {input.non_host_r1} \
                -2 {input.non_host_r2} \
                --input-profile {output.condense} \
                --output-tsv {output.read_fraction} \
                --output-per-taxon-read-fractions {params.read_fraction_taxa}
            fi

        #Otheriwse, don't run singlem
        else
        echo "SingleM analysis not performed"
        touch {output.condense}
        touch {output.pipe}
        echo -e "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\n0\t0\t0\t0.0%" > {output.read_fraction}
        
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