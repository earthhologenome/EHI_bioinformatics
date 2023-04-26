################################################################################
################################################################################
################################################################################
# EHI snakefile for preprocessing raw reads (trimming/mapping to host)
# Raphael Eisenhofer 02/2023
#         .----------------.  .----------------.  .----------------.
#        | .--------------. || .--------------. || .--------------. |
#        | |  _________   | || |  ____  ____  | || |     _____    | |
#        | | |_   ___  |  | || | |_   ||   _| | || |    |_   _|   | |
#        | |   | |_  \_|  | || |   | |__| |   | || |      | |     | |
#        | |   |  _|  _   | || |   |  __  |   | || |      | |     | |
#        | |  _| |___/ |  | || |  _| |  | |_  | || |     _| |_    | |
#        | | |_________|  | || | |____||____| | || |    |_____|   | |
#        | |              | || |              | || |              | |
#        | '--------------' || '--------------' || '--------------' |
#         '----------------'  '----------------'  '----------------'
################################################################################
################################################################################
################################################################################

configfile: "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/configs/1_Preprocess_QC_config.yaml"

### Setup sample inputs, config, and working directory
import pandas as pd

SAMPLE = pd.read_csv('PRBATCH_input.tsv', sep='\t', header=None).loc[:, 0].tolist()

print("Detected the following samples:")
print(SAMPLE)


### Code to scale time needed by raw read file sizes
### Scaling is based on benchmark data for ~280 jobs 3/4/2023 RE
import os

def estimate_time_download(wildcards):
    fs_sample = f"/projects/ehi/data/RAW/PRBATCH/{wildcards.sample}_filesize.txt"
    with open(fs_sample, 'r') as f:
        input_size = int(f.read().strip())
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # Multiply by 2, and set time based on 30 MB/s download speed.
    estimate_time_download = ((input_size_gb * 2.1 ) + 12 ) / 1.25
    return int(estimate_time_download)

def estimate_time_fastp(wildcards):
    r1_path = f"/projects/ehi/data/RAW/PRBATCH/{wildcards.sample}_1.fq.gz"
    r2_path = f"/projects/ehi/data/RAW/PRBATCH/{wildcards.sample}_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # Add scaling (* 2.5 is for the Gbp to .gz compressed filesize scaling -- e.g. 3 Gbp sample ~ 1.5 GBytes) 
    estimate_time_fastp = ((input_size_gb * 2.5 ) + 6) / 2
    return int(estimate_time_fastp)

def estimate_time_mapping(wildcards):
    r1_path = f"/projects/ehi/data/PPR/PRBATCH/tmp/{wildcards.sample}_trimmed_1.fq.gz"
    r2_path = f"/projects/ehi/data/PPR/PRBATCH/tmp/{wildcards.sample}_trimmed_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    estimate_time_mapping = ((input_size_gb * 2 ) + 2) * 14
    return int(estimate_time_mapping)

def estimate_time_nonpareil(wildcards):
    r1_path = f"/projects/ehi/data/PPR/PRBATCH/{wildcards.sample}_M_1.fq"
    r2_path = f"/projects/ehi/data/PPR/PRBATCH/{wildcards.sample}_M_2.fq"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    # N.b. for estimate_time_nonpareil we estimate from uncompressed fq
    estimate_time_nonpareil = (input_size_gb + 2) * 2
    return int(estimate_time_nonpareil)


def estimate_time_singlem(wildcards):
    r1_path = f"/projects/ehi/data/PPR/PRBATCH/{wildcards.sample}_M_1.fq.gz"
    r2_path = f"/projects/ehi/data/PPR/PRBATCH/{wildcards.sample}_M_2.fq.gz"
    input_files = [r1_path, r2_path]
    input_size = sum(os.path.getsize(f) for f in input_files)
    # convert from bytes to gigabytes
    input_size_gb = input_size / (1024 * 1024 * 1024)
    estimate_time_singlem = ((input_size_gb) + 7) * 5
    return int(estimate_time_singlem)

################################################################################
### Setup the desired outputs
rule all:
    input:
        "/projects/ehi/data/REP/PRBATCH.tsv",
        "/projects/ehi/data/PPR/PRBATCH/0_REPORTS/PRBATCH_nonpareil_metadata.tsv",
################################################################################
### Create PRB folder on ERDA
rule create_PRB_folder:
    output:
        "/projects/ehi/data/PPR/PRBATCH/ERDA_folder_created"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=1,
        mem_gb=8,
        time='00:03:00'
    message:
        "Creating PRB folder on ERDA"
    shell:
        """
        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/PPR/PRBATCH ; bye"
        touch {output}

        #Also, log the AirTable that the PRB is running!
        python /projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/airtable/log_prb_start_airtable.py --code=PRBATCH

        """
################################################################################
### Get file sizes from ERDA
rule filesize_from_ERDA:
    input:
        "/projects/ehi/data/PPR/PRBATCH/ERDA_folder_created"
    output:
        temp("/projects/ehi/data/RAW/PRBATCH/{sample}_filesize.txt"),
    params:
        workdir = config["workdir"]
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time='00:00:30'
    message:
        "Fetching filesize for {wildcards.sample} from ERDA"
    shell:
        """
        echo 'ls -l /EarthHologenomeInitiative/Data/RAW/*/{wildcards.sample}*_1.fq.gz' | sftp erda | sed '1d;' | awk '{{print $5}}' > {output}

        """
################################################################################
### Fetch raw data from ERDA
rule download_from_ERDA:
    input:
        ready = "/projects/ehi/data/PPR/PRBATCH/ERDA_folder_created",
        filesize = "/projects/ehi/data/RAW/PRBATCH/{sample}_filesize.txt"
    output:
        r1o = temp("/projects/ehi/data/RAW/PRBATCH/{sample}_1.fq.gz"),
        r2o = temp("/projects/ehi/data/RAW/PRBATCH/{sample}_2.fq.gz"),
    params:
        workdir = config["workdir"]
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time=estimate_time_download
    message:
        "Fetching {wildcards.sample} from ERDA"
    shell:
        """
        lftp sftp://erda -e "mirror --include-glob='{wildcards.sample}*.fq.gz' /EarthHologenomeInitiative/Data/RAW/ {params.workdir}/RAW/PRBATCH/; bye"
        mv {params.workdir}/RAW/PRBATCH/*/{wildcards.sample}*_1.fq.gz {output.r1o}
        mv {params.workdir}/RAW/PRBATCH/*/{wildcards.sample}*_2.fq.gz {output.r2o}
        """
################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i = "/projects/ehi/data/RAW/PRBATCH/{sample}_1.fq.gz",
        r2i = "/projects/ehi/data/RAW/PRBATCH/{sample}_2.fq.gz",
    output:
        r1o = temp("/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_1.fq.gz"),
        r2o = temp("/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_2.fq.gz"),
        fastp_html = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}.html",
        fastp_json = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}.json"
    params:
        adapter1 = expand("{adapter1}", adapter1=config['adapter1']),
        adapter2 = expand("{adapter2}", adapter2=config['adapter2'])
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=24,
        time=estimate_time_fastp
    benchmark:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_fastp.benchmark.tsv"
    log:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
            --low_complexity_filter \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence {params.adapter1} \
            --adapter_sequence_r2 {params.adapter2} \
        &> {log}
        """
################################################################################
## Fetch host genome from ERDA, if not there already, download and index it.
rule fetch_host_genome:
    input:
        "/projects/ehi/data/PPR/PRBATCH/ERDA_folder_created"
    output:
        bt2_index = "/projects/ehi/data/RUN/PRBATCH/HOST_GENOME/HOST_GENOME_RN.fna.gz.rev.2.bt2l",
        rn_catted_ref = "/projects/ehi/data/RUN/PRBATCH/HOST_GENOME/HOST_GENOME_RN.fna.gz"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    params:
        workdir = config["workdir"]
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='03:00:00'
    log:
        "/projects/ehi/data/RUN/PRBATCH/logs/host_genome_indexing.log"
    message:
        "Fetching host genome"
    shell:
        """
        # IF statement for if file exists on Mjolnir
        if [ -f {output.bt2_index} ]
            then
                echo "Genome is ready to go!"

            elif 
                sftp_check=$(sftp erda:/EarthHologenomeInitiative/Data/GEN/HOST_GENOME.tar.gz 2>&1)
                echo "$sftp_check" | grep -q "not found"

            then
                echo "Downloading and indexing reference genome"
                mkdir -p {params.workdir}/RUN/PRBATCH/HOST_GENOME/
                wget HG_URL -q -O {params.workdir}/RUN/PRBATCH/HOST_GENOME/HOST_GENOME.fna.gz

                # Add '_' separator for CoverM
                rename.sh \
                    in={params.workdir}/RUN/PRBATCH/HOST_GENOME/HOST_GENOME.fna.gz \
                    out={output.rn_catted_ref} \
                    prefix=HOST_GENOME \
                    -Xmx{resources.mem_gb}G 
                
                rm {params.workdir}/RUN/PRBATCH/HOST_GENOME/HOST_GENOME.fna.gz

                # Index catted genomes
                bowtie2-build \
                    --large-index \
                    --threads {threads} \
                    {output.rn_catted_ref} {output.rn_catted_ref} \
                    &> {log}

                # Compress and upload to ERDA for future use
                cd {params.workdir}/RUN/PRBATCH/HOST_GENOME/
                tar -I pigz -cvf HOST_GENOME.tar.gz *
                sftp erda:/EarthHologenomeInitiative/Data/GEN/ <<< $'put HOST_GENOME.tar.gz'
                rm HOST_GENOME.tar.gz
                cd {params.workdir}/RUN/PRBATCH

                # Log AirTable that a new genome has been indexed and uploaded to ERDA
                python /projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/airtable/log_genome_airtable.py --code=HOST_GENOME

            else 
                echo "Indexed genome exists on erda, unpacking."
                tar -xvzf HOST_GENOME.tar.gz --directory {params.workdir}/RUN/PRBATCH/HOST_GENOME/
                rm HOST_GENOME.tar.gz

        fi

        """

################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i = "/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_1.fq.gz",
        r2i = "/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_2.fq.gz",
        catted_ref = "/projects/ehi/data/RUN/PRBATCH/HOST_GENOME/HOST_GENOME_RN.fna.gz",
        bt2_index = "/projects/ehi/data/RUN/PRBATCH/HOST_GENOME/HOST_GENOME_RN.fna.gz.rev.2.bt2l"
    output:
        all_bam = temp("/projects/ehi/data/PPR/PRBATCH/tmp/{sample}.bam"),
        host_bam = temp("/projects/ehi/data/PPR/PRBATCH/{sample}_G.bam"),
        non_host_r1 = temp("/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq"),
        non_host_r2 = temp("/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq"),
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=24,
        time=estimate_time_mapping
    benchmark:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_mapping.benchmark.tsv"
    log:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_mapping.log"
    message:
        "Mapping {wildcards.sample} reads to host genomes"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.catted_ref} \
            -1 {input.r1i} \
            -2 {input.r2i} \
        | samtools view -b -@ {threads} - | samtools sort -@ {threads} -o {output.all_bam} - &&

        # Extract non-host reads (note we're not compressing for nonpareil)
        samtools view -b -f12 -@ {threads} {output.all_bam} \
        | samtools fastq -@ {threads} -1 {output.non_host_r1} -2 {output.non_host_r2} - &&

        # Send host reads to BAM
        samtools view -b -F12 -@ {threads} {output.all_bam} \
        | samtools sort -@ {threads} -o {output.host_bam} -
        """




################################################################################
### Estimate diversity and required sequencing effort using nonpareil
rule nonpareil:
    input:
        non_host_r1 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq",
        non_host_r2 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq",
    output:
        npo = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}.npo",
        non_host_r1c = temp("/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq.gz"),
        non_host_r2c = temp("/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq.gz"),
    params:
        sample = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}",
        workdir = config["workdir"]
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/nonpareil.yaml"
    threads:
        8
    resources:
        load=1,
        mem_gb=45,
        time=estimate_time_nonpareil
    benchmark:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_nonpareil.benchmark.tsv"
    message:
        "Estimating microbial diversity using nonpareil"
    shell:
        """
        #IF statement to account for situations where there are not enough
        #microbial reads in a sample (e.g. high host% or non-metagenomic sample)
        #In this case, if R1 has > 150 Mbytes, run, else, skip:
        if [ $(( $(stat -c '%s' {input.non_host_r1}) / 1024 / 1024 )) -gt 150 ]
        then
        #Run nonpareil
        nonpareil \
            -s {input.non_host_r1} \
            -f fastq \
            -T kmer \
            -t {threads} \
            -b {params.sample}
        else
        #Create dummy file for snakemake to proceed
        touch {output.npo}
        fi

        #Compress reads
        pigz -p {threads} {input.non_host_r1}
        pigz -p {threads} {input.non_host_r2}

        """
################################################################################
### Estimate the fraction of bacterial and archaeal DNA using SingleM read fraction
rule singlem:
    input:
        npo = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}.npo",
        non_host_r1 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq.gz",
    output:
        pipe = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_pipe.tsv.gz",
        condense = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_condense.tsv",
        read_fraction = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_readfraction.tsv",
    params:
        pipe_uncompressed = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_pipe.tsv",
        read_fraction_taxa = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_readfraction_per_taxa.tsv"
# Current issue with snakemake and pre-built conda environments: https://github.com/snakemake/snakemake/pull/1708
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/singlem.yaml"
    threads:
        3
    resources:
        load=1,
        mem_gb=36,
        time=estimate_time_singlem
    benchmark:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_singlem.benchmark.tsv"
    message:
        "Estimating microbial fraction using singlem"
    shell:
        """
        #Temp fix until snakemake is fixed or singlem conda recipe is updated
        export PATH='/projects/ehi/data/0_Environments/github_repos/singlem/bin':$PATH
        export SINGLEM_METAPACKAGE_PATH='/projects/ehi/data/0_Environments/databases/S3.1.0.metapackage_20221209.smpkg.zb/'

        #Try to fix /tmp folder running out of space:
        mkdir -p TMPDIR=/projects/ehi/data/PPR/PRBATCH/tmp/tmp
        export TMPDIR=/projects/ehi/data/PPR/PRBATCH/tmp/tmp

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
################################################################################
### Calculate % of each sample's reads mapping to host genome/s (also upload PPR reads to ERDA)
rule coverM:
    input:
        bam = "/projects/ehi/data/PPR/PRBATCH/{sample}_G.bam",
        npo = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}.npo",
        pipe = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_pipe.tsv.gz",
        non_host_r1 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq.gz",
        host_bam = "/projects/ehi/data/PPR/PRBATCH/{sample}_G.bam",
    output:
        "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_coverM_mapped_host.tsv"
    params:
        assembly = "/projects/ehi/data/RUN/PRBATCH/HOST_GENOME/HOST_GENOME_RN.fna.gz",
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/coverm.yaml"
    threads:
        2
    resources:
        load=1,
        mem_gb=16,
        time='00:10:00'
    benchmark:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_coverM.benchmark.tsv"
    log:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_coverM.log"
    message:
        "Calculating percentage of reads mapped to host genome/s using coverM"
    shell:
        """
        #Calculate % mapping to host using coverM
        coverm genome \
            -b {input.bam} \
            -s _ \
            -m relative_abundance count \
            -t {threads} \
            --min-covered-fraction 0 \
            > {output}

        #Remove empty nonpareil (.npo) files to streamline future plotting
        if [ $(stat -c '%s' {input.npo}) -lt 1 ]
        then
        rm {input.npo}
        fi

        """
################################################################################
### Calculate % of each sample's reads mapping to host genome/s (also upload PPR reads to ERDA)
rule upload_to_ERDA:
    input:
        non_host_r1 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2 = "/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq.gz",
        host_bam = "/projects/ehi/data/PPR/PRBATCH/{sample}_G.bam",
        coverm = "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_coverM_mapped_host.tsv",
        fiel_size = "/projects/ehi/data/RAW/PRBATCH/{sample}_filesize.txt",
    output:
        "/projects/ehi/data/PPR/PRBATCH/misc/{sample}_uploaded"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=16,
        time=estimate_time_download
    log:
        "/projects/ehi/data/RUN/PRBATCH/logs/{sample}_upload.log"
    message:
        "Uploading reads and BAM to ERDA"
    shell:
        """
        #Upload preprocessed reads to ERDA for storage
        lftp sftp://erda -e "put {input.non_host_r1} -o /EarthHologenomeInitiative/Data/PPR/PRBATCH/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.non_host_r2} -o /EarthHologenomeInitiative/Data/PPR/PRBATCH/; bye"
        sleep 5
        lftp sftp://erda -e "put {input.host_bam} -o /EarthHologenomeInitiative/Data/PPR/PRBATCH/; bye"
        touch {output}
        """
################################################################################
### Create summary table from outputs
rule report:
    input:
        coverm = expand("/projects/ehi/data/PPR/PRBATCH/misc/{sample}_coverM_mapped_host.tsv", sample=SAMPLE),
        fastp = expand("/projects/ehi/data/PPR/PRBATCH/misc/{sample}.json", sample=SAMPLE),
        read_fraction = expand("/projects/ehi/data/PPR/PRBATCH/misc/{sample}_readfraction.tsv", sample=SAMPLE),
        file_size = expand("/projects/ehi/data/RAW/PRBATCH/{sample}_filesize.txt", sample=SAMPLE),
        uploaded = expand("/projects/ehi/data/PPR/PRBATCH/misc/{sample}_uploaded", sample=SAMPLE),
        raw1 = expand("/projects/ehi/data/RAW/PRBATCH/{sample}_1.fq.gz", sample=SAMPLE),
        raw2 = expand("/projects/ehi/data/RAW/PRBATCH/{sample}_2.fq.gz", sample=SAMPLE),
        trimmed1 = expand("/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_1.fq.gz", sample=SAMPLE),
        trimmed2 = expand("/projects/ehi/data/PPR/PRBATCH/tmp/{sample}_trimmed_2.fq.gz", sample=SAMPLE),
        non_host_r1 = expand("/projects/ehi/data/PPR/PRBATCH/{sample}_M_1.fq.gz", sample=SAMPLE),
        non_host_r2 = expand("/projects/ehi/data/PPR/PRBATCH/{sample}_M_2.fq.gz", sample=SAMPLE),
    output:
        report = "/projects/ehi/data/REP/PRBATCH.tsv",
        npar_metadata = "/projects/ehi/data/PPR/PRBATCH/0_REPORTS/PRBATCH_nonpareil_metadata.tsv"
    params:
        tmpdir = "/projects/ehi/data/PPR/PRBATCH/tmp/",
        npar = expand("/projects/ehi/data/PPR/PRBATCH/misc/{sample}.npo", sample=SAMPLE),
        misc_dir = "/projects/ehi/data/PPR/PRBATCH/misc/",
        workdir = config["workdir"]
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=1,
        mem_gb=24,
        time='00:20:00'
    message:
        "Creating a final preprocessing report"
    shell:
        """
        #Create nonpareil sample metadata file
        mkdir -p {params.tmpdir}
        for i in {params.npar}; do echo $(basename $i) >> {params.tmpdir}/files.txt; done
        for i in {params.npar}; do echo $(basename ${{i/.npo/}}) >> {params.tmpdir}/names.txt; done
        for i in {params.npar}; do echo "#f03b20" >> {params.tmpdir}/colours.txt; done
        echo -e "File\tName\tColour" > {params.tmpdir}/headers.txt
        paste {params.tmpdir}/files.txt {params.tmpdir}/names.txt {params.tmpdir}/colours.txt > {params.tmpdir}/merged.tsv
        cat {params.tmpdir}/headers.txt {params.tmpdir}/merged.tsv > {output.npar_metadata}

        #Create preprocessing report
        mkdir -p {params.tmpdir}
        for i in {input.coverm}; do echo $(basename ${{i/_coverM_mapped_host.tsv}}) >> {params.tmpdir}/names.tsv; done
        for i in {input.coverm}; do grep -v 'Genome' $i | grep -v 'unmapped' | cut -f3; done >> {params.tmpdir}/host_reads.tsv

        for i in {input.fastp}; do grep '"total_reads"' $i | sed -n 1p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/read_pre_filt.tsv
        for i in {input.fastp}; do grep '"total_reads"' $i | sed -n 2p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/read_post_filt.tsv
        for i in {input.fastp}; do grep '"total_bases"' $i | sed -n 1p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/bases_pre_filt.tsv
        for i in {input.fastp}; do grep '"total_bases"' $i | sed -n 2p | cut -f2 --delimiter=: | tr -d ','; done >> {params.tmpdir}/bases_post_filt.tsv
        for i in {input.fastp}; do grep 'adapter_trimmed_reads' $i | cut -f2 --delimiter=: | tr -d ',' | tr -d ' '; done >> {params.tmpdir}/adapter_trimmed_reads.tsv
        for i in {input.fastp}; do grep 'adapter_trimmed_bases' $i | cut -f2 --delimiter=: | tr -d ',' | tr -d ' '; done >> {params.tmpdir}/adapter_trimmed_bases.tsv

        #parse singlem estimates
        for i in {input.read_fraction}; do sed '1d;' $i | cut -f2,3,4 >> {params.tmpdir}/singlem.tsv; done

        paste {params.tmpdir}/names.tsv {params.tmpdir}/read_pre_filt.tsv {params.tmpdir}/read_post_filt.tsv {params.tmpdir}/bases_pre_filt.tsv {params.tmpdir}/bases_post_filt.tsv {params.tmpdir}/adapter_trimmed_reads.tsv {params.tmpdir}/adapter_trimmed_bases.tsv {params.tmpdir}/host_reads.tsv {params.tmpdir}/singlem.tsv > {params.tmpdir}/preprocessing_stats.tsv
        echo -e "EHI_number\treads_pre_fastp\treads_post_fastp\tbases_pre_fastp\tbases_post_fastp\tadapter_trimmed_reads\tadapter_trimmed_bases\thost_reads\tbacterial_archaeal_bases\tmetagenomic_bases\tsinglem_fraction" > {params.tmpdir}/headers.tsv
        cat {params.tmpdir}/headers.tsv {params.tmpdir}/preprocessing_stats.tsv > {output.report}

        cp {output.report} {params.misc_dir}
        cp {output.npar_metadata} {params.misc_dir}
        tar -czf PRBATCH_stats.tar.gz {params.misc_dir}

        #Upload stats and report to ERDA for storage
        lftp sftp://erda -e "put PRBATCH_stats.tar.gz -o /EarthHologenomeInitiative/Data/PPR/PRBATCH/; bye"
        sleep 10
        lftp sftp://erda -e "put {output.report} -o /EarthHologenomeInitiative/Data/REP/; bye"

        #Automatically update the AirTable with the preprocessing stats
        python /projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/airtable/add_prb_stats_airtable.py --report={output.report} --prb=PRBATCH 

        #Indicate that the PRB is done in AirTable
        python /projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/airtable/log_prb_done_airtable.py --code=PRBATCH
       
        #Clean up the files/directories
        rm PRBATCH_stats.tar.gz
        rm -r {params.workdir}/RUN/PRBATCH/HOST_GENOME/
        rm -r {params.misc_dir}
        rm -r {params.tmpdir}

        """
