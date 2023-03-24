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

################################################################################
### Setup the desired outputs
rule all:
    input:
        "$workdir/REP/PRBATCH.tsv",
        "$workdir/PPR/PRBATCH/0_REPORTS/PRBATCH_nonpareil_metadata.tsv",
        "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz.rev.2.bt2l"
################################################################################
### Fetch raw data from ERDA
rule fetch_raw_reads:
    output:
        r1o = temp("$workdir/RAW/PRBATCH/{sample}_1.fq.gz"),
        r2o = temp("$workdir/RAW/PRBATCH/{sample}_2.fq.gz"),
    params:
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        1
    resources:
        mem_gb=8,
        time='00:15:00'
    message:
        "Fetching {wildcards.sample} from ERDA"
    shell:
        """
        sftp erda:/EarthHologenomeInitiative/Data/RAW/*/{wildcards.sample}*.fq.gz $workdir/RAW/PRBATCH/
        mv $workdir/RAW/PRBATCH/{wildcards.sample}*_1.fq.gz {output.r1o}
        mv $workdir/RAW/PRBATCH/{wildcards.sample}*_2.fq.gz {output.r2o}
        """
################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i = "$workdir/RAW/PRBATCH/{sample}_1.fq.gz",
        r2i = "$workdir/RAW/PRBATCH/{sample}_2.fq.gz"
    output:
        r1o = temp("$workdir/PPR/PRBATCH/tmp/{sample}_trimmed_1.fq.gz"),
        r2o = temp("$workdir/PPR/PRBATCH/tmp/{sample}_trimmed_2.fq.gz"),
        fastp_html = "$workdir/PPR/PRBATCH/misc/{sample}.html",
        fastp_json = "$workdir/PPR/PRBATCH/misc/{sample}.json"
    params:
        adapter1 = expand("{adapter1}", adapter1=config['adapter1']),
        adapter2 = expand("{adapter2}", adapter2=config['adapter2'])
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        mem_gb=24,
        time='00:30:00'
    benchmark:
        "$workdir/RUN/PRBATCH/logs/{sample}_fastp.benchmark.tsv"
    log:
        "$workdir/RUN/PRBATCH/logs/{sample}_fastp.log"
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
    output:
        bt2_index = "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz.rev.2.bt2l",
        rn_catted_ref = "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        16
    resources:
        mem_gb=96,
        time='03:00:00'
    log:
        "$workdir/PRBATCH/logs/host_genome_indexing.log"
    message:
        "Fetching host genome"
    shell:
        """
        # IF statement for if file exists on Mjolnir
        if [ -f {output.rn_catted_ref} ]
            then
                echo "Genome is ready to go!"

            elif 
                sftp_check=$(sftp erda:/EarthHologenomeInitiative/Data/GEN/HOST_GENOME.tar.gz 2>&1)
                echo "$sftp_check" | grep -q "not found"

            then
                echo "Downloading and indexing reference genome"
                mkdir -p $workdir/GEN/HOST_GENOME/
                wget HG_URL -q -O $workdir/GEN/HOST_GENOME/HOST_GENOME.fna.gz

                # Add '_' separator for CoverM
                rename.sh \
                    in=$workdir/GEN/HOST_GENOME/HOST_GENOME.fna.gz \
                    out={output.rn_catted_ref} \
                    prefix=HOST_GENOME \
                    -Xmx{resources.mem_gb}G 
                
                rm $workdir/GEN/HOST_GENOME/HOST_GENOME.fna.gz

                # Index catted genomes
                bowtie2-build \
                    --large-index \
                    --threads {threads} \
                    {output.rn_catted_ref} {output.rn_catted_ref} \
                    &> {log}

                # Compress and upload to ERDA for future use
                tar -I pigz -cvf $workdir/GEN/HOST_GENOME/HOST_GENOME.tar.gz $workdir/GEN/HOST_GENOME/*
                sftp erda:/EarthHologenomeInitiative/Data/GEN/ <<< $'put $workdir/GEN/HOST_GENOME/HOST_GENOME.tar.gz'
                rm $workdir/GEN/HOST_GENOME/HOST_GENOME.tar.gz

                # Create a warning that a new genome has been indexed and needs to be logged in AirTable
                echo "HOST_GENOME has been indexed and needs to be logged in AirTable" > $workdir/NEW_HOST_GENOME.txt
                
            else 
                echo "Indexed genome exists on erda, unpacking."
                mv HOST_GENOME.tar.gz $workdir/
                tar -xvzf $workdir/HOST_GENOME.tar.gz
                rm $workdir/HOST_GENOME.tar.gz

        fi

        # Create PRBATCH folder on ERDA for uploading processed reads and BAMs
        sftp erda:/EarthHologenomeInitiative/Data/PPR <<< $'mkdir PRBATCH'

        # tmp file for if another person is using the same genome? so it doesn't get deleted?
        """
################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i = "$workdir/PPR/PRBATCH/tmp/{sample}_trimmed_1.fq.gz",
        r2i = "$workdir/PPR/PRBATCH/tmp/{sample}_trimmed_2.fq.gz",
        catted_ref = "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz",
        bt2_index = "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz.rev.2.bt2l"
    output:
        all_bam = temp("$workdir/PPR/PRBATCH/tmp/{sample}.bam"),
        host_bam = temp("$workdir/PPR/PRBATCH/{sample}_G.bam"),
        non_host_r1 = temp("$workdir/PPR/PRBATCH/{sample}_M_1.fq"),
        non_host_r2 = temp("$workdir/PPR/PRBATCH/{sample}_M_2.fq"),
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        mem_gb=48,
        time='03:00:00'
    benchmark:
        "$workdir/RUN/PRBATCH/logs/{sample}_mapping.benchmark.tsv"
    log:
        "$workdir/RUN/PRBATCH/logs/{sample}_mapping.log"
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
        non_host_r1 = "$workdir/PPR/PRBATCH/{sample}_M_1.fq",
        non_host_r2 = "$workdir/PPR/PRBATCH/{sample}_M_2.fq",
    output:
        npo = "$workdir/PPR/PRBATCH/misc/{sample}.npo",
        non_host_r1c = "$workdir/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2c = "$workdir/PPR/PRBATCH/{sample}_M_2.fq.gz",
    params:
        sample = "$workdir/PPR/PRBATCH/misc/{sample}",
        badsample_r1 = "$workdir/PPR/PRBATCH/poor_samples/{sample}_M_1.fq.gz",
        badsample_r2 = "$workdir/PPR/PRBATCH/poor_samples/{sample}_M_2.fq.gz"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/nonpareil.yaml"
    threads:
        8
    resources:
        mem_gb=45,
        time='02:00:00'
    benchmark:
        "$workdir/RUN/PRBATCH/logs/{sample}_nonpareil.benchmark.tsv"
    message:
        "Estimating microbial diversity using nonpareil"
    shell:
        """
        mkdir -p $workdir/PPR/PRBATCH/poor_samples

        #IF statement to account for situations where there are not enough
        #microbial reads in a sample (e.g. high host% or non-metagenomic sample)
        #In this case, if R1 has > 100 Mbytes, run, else, skip:
        if [ $(( $(stat -c '%s' {input.non_host_r1}) / 1024 / 1024 )) -gt 100 ]
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

        #Move samples that don't have enough reads for assembly to a new folder
        #This saves time, and prevents errors in the next pipeline!
        if [ $(( $(stat -c '%s' {input.non_host_r1}.gz) / 1024 / 1024 )) -lt 200 ]
        then
        mv {input.non_host_r1}.gz {params.badsample_r1} && mv {input.non_host_r2}.gz {params.badsample_r2}
        fi
        """
################################################################################
### Estimate the fraction of bacterial and archaeal DNA using SingleM read fraction
rule singlem:
    input:
        non_host_r1 = "$workdir/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2 = "$workdir/PPR/PRBATCH/{sample}_M_2.fq.gz",
    output:
        pipe = "$workdir/PPR/PRBATCH/misc/{sample}_pipe.tsv.gz",
        condense = "$workdir/PPR/PRBATCH/misc/{sample}_condense.tsv",
        read_fraction = "$workdir/PPR/PRBATCH/misc/{sample}_readfraction.tsv",
    params:
        pipe_uncompressed = "$workdir/PPR/PRBATCH/misc/{sample}_pipe.tsv",
        read_fraction_taxa = "$workdir/PPR/PRBATCH/misc/{sample}_readfraction_per_taxa.tsv"
# Current issue with snakemake and pre-built conda environments: https://github.com/snakemake/snakemake/pull/1708
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/singlem.yaml"
    threads:
        8
    resources:
        mem_gb=45,
        time='02:00:00'
    benchmark:
        "$workdir/RUN/PRBATCH/logs/{sample}_singlem.benchmark.tsv"
    message:
        "Estimating microbial fraction using singlem"
    shell:
        """
        #Temp fix until snakemake is fixed or singlem conda recipe is updated
        export PATH='/projects/ehi/data/0_Environments/github_repos/singlem/bin':$PATH
        
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
        echo -e "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\nNA\tNA\tNA\tNA" > {output.read_fraction}
        
        else        
        #Run singlem read_fraction
        singlem read_fraction \
            -1 {input.non_host_r1} \
            -2 {input.non_host_r2} \
            --input-profile {output.condense} \
            --output-tsv {output.read_fraction} \
            --output-per-taxon-read-fractions {params.read_fraction_taxa}
        fi
        
        #If statement for cases when singlem does not produce a condense output
        if [ -f {params.read_fraction_taxa} ]
        then
        #Compress read_fraction_per_taxa file
        gzip {params.read_fraction_taxa}

        else
        echo "no microbes in sample"
        fi
        
        """
################################################################################
### Calculate % of each sample's reads mapping to host genome/s (also upload PPR reads to ERDA)
rule coverM_and_upload_to_ERDA:
    input:
        bam = "$workdir/PPR/PRBATCH/{sample}_G.bam",
        npo = "$workdir/PPR/PRBATCH/misc/{sample}.npo",
        pipe = "$workdir/PPR/PRBATCH/misc/{sample}_pipe.tsv.gz",
        non_host_r1 = "$workdir/PPR/PRBATCH/{sample}_M_1.fq.gz",
        non_host_r2 = "$workdir/PPR/PRBATCH/{sample}_M_2.fq.gz",
        host_bam = "$workdir/PPR/PRBATCH/{sample}_G.bam",
    output:
        "$workdir/PPR/PRBATCH/misc/{sample}_coverM_mapped_host.tsv"
    params:
        assembly = "$workdir/GEN/HOST_GENOME/HOST_GENOME_RN.fna.gz",
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/coverm.yaml"
    threads:
        2
    resources:
        mem_gb=16,
        time='01:00:00'
    benchmark:
        "$workdir/RUN/PRBATCH/logs/{sample}_coverM.benchmark.tsv"
    log:
        "$workdir/RUN/PRBATCH/logs/{sample}_coverM.log"
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

        #Upload preprocessed reads to ERDA for storage
        sftp erda:/EarthHologenomeInitiative/Data/PPR/PRBATCH/ <<< $'put {input.non_host_r1}'
        sftp erda:/EarthHologenomeInitiative/Data/PPR/PRBATCH/ <<< $'put {input.non_host_r2}'
        sftp erda:/EarthHologenomeInitiative/Data/PPR/PRBATCH/ <<< $'put {input.host_bam}'
        """
################################################################################
### Create summary table from outputs
rule report:
    input:
        coverm = expand("$workdir/PPR/PRBATCH/misc/{sample}_coverM_mapped_host.tsv", sample=SAMPLE),
        fastp = expand("$workdir/PPR/PRBATCH/misc/{sample}.json", sample=SAMPLE),
        read_fraction = expand("$workdir/PPR/PRBATCH/misc/{sample}_readfraction.tsv", sample=SAMPLE)
    output:
        report = "$workdir/REP/PRBATCH.tsv",
        npar_metadata = "$workdir/PPR/PRBATCH/0_REPORTS/PRBATCH_nonpareil_metadata.tsv"
    params:
        tmpdir = "$workdir/PPR/PRBATCH/tmp/",
        npar = expand("$workdir/PPR/PRBATCH/misc/{sample}.npo", sample=SAMPLE),
        misc_dir = "$workdir/PPR/PRBATCH/misc/"
    conda:
        "/projects/ehi/data/0_Code/EHI_bioinformatics_EHI_VERSION/0_Code/conda_envs/1_Preprocess_QC.yaml"
    threads:
        1
    resources:
        mem_gb=45,
        time='00:05:00'
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
        echo -e "sample\treads_pre_filt\treads_post_filt\tbases_pre_filt\tbases_post_filt\tadapter_trimmed_reads\tadapter_trimmed_bases\thost_reads\tbacterial_archaeal_bases\tmetagenomic_bases\tsinglem_fraction" > {params.tmpdir}/headers.tsv
        cat {params.tmpdir}/headers.tsv {params.tmpdir}/preprocessing_stats.tsv > {output.report}

        mv slurm*

        cp {output.report} {params.misc_dir}
        cp {output.npar_metadata} {params.misc_dir}
        tar -czf PRBATCH_stats.tar.gz {params.misc_dir}

        rm -r {params.tmpdir}

        #Upload stats and report to ERDA for storage
        sftp erda:/EarthHologenomeInitiative/Data/PPR/PRBATCH/ <<< $'put PRBATCH_stats.tar.gz'
        sftp erda:/EarthHologenomeInitiative/Data/RUN/ <<< $'put {output.report}'

        #Clean up the files/directories
        rm PRBATCH_stats.tar.gz
        rm -r $workdir/GEN/HOST_GENOME/

        """
