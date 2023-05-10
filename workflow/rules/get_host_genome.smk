################################################################################
## Fetch host genome from ERDA, if not there already, download and index it.
rule fetch_host_genome:
    input:
        os.path.join(
            config["workdir"], 
            "ERDA_folder_created"
        )
    output:
        bt2_index=os.path.join(
            config["workdir"],
            config["hostgenome"],
            config["hostgenome"] + "_RN.fna.gz.rev.2.bt2l",
        ),
        rn_catted_ref=os.path.join(
            config["workdir"],
            config["hostgenome"],
            config["hostgenome"] + "_RN.fna.gz"
        )
    conda:
        f"{config['codedir']}/conda_envs/1_Preprocess_QC.yaml"
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='03:00:00'
    log:
        os.path.join(config["logdir"] + "/host_genome_indexing.log")
    message:
        "Fetching host genome"
    shell:
        """
        # IF statement for if file exists on Mjolnir
        if [ -f {output.bt2_index} ]
            then
                echo "Genome is ready to go!"

            elif 
                sftp_check=$(sftp erda:/EarthHologenomeInitiative/Data/GEN/{config[hostgenome]}.tar.gz 2>&1)
                echo "$sftp_check" | grep -q "not found"

            then
                echo "Downloading and indexing reference genome"
                mkdir -p {config[workdir]}/{config[hostgenome]}/
                wget {config[hg_url]} -q -O {config[workdir]}/{config[hostgenome]}/{config[hostgenome]}.fna.gz

                # Add '_' separator for CoverM
                rename.sh \
                    in={config[workdir]}/{config[hostgenome]}/{config[hostgenome]}.fna.gz \
                    out={output.rn_catted_ref} \
                    prefix={config[hostgenome]} \
                    -Xmx{resources.mem_gb}G 
                
                rm {config[workdir]}/{config[hostgenome]}/{config[hostgenome]}.fna.gz

                # Index catted genomes
                bowtie2-build \
                    --large-index \
                    --threads {threads} \
                    {output.rn_catted_ref} {output.rn_catted_ref} \
                    &> {log}

                # Compress and upload to ERDA for future use
                cd {config[workdir]}/{config[hostgenome]}/
                tar -I pigz -cvf {config[hostgenome]}.tar.gz *
                sftp erda:/EarthHologenomeInitiative/Data/GEN/ <<< $'put {config[hostgenome]}.tar.gz'
                rm {config[hostgenome]}.tar.gz
                cd {config[workdir]}

                # Log AirTable that a new genome has been indexed and uploaded to ERDA
                python {config[codedir]}/airtable/log_genome_airtable.py --code={config[hostgenome]}

            else 
                echo "Indexed genome exists on erda, unpacking."
                tar -xvzf {config[hostgenome]}.tar.gz --directory {config[workdir]}/{config[hostgenome]}/
                rm {config[hostgenome]}.tar.gz

        fi

        """