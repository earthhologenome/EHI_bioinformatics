################################################################################
################################################################################
################################################################################
# EHI snakefile for preprocessing raw reads (trimming/mapping to host)
# Raphael Eisenhofer 4/2022
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

configfile: "0_Code/1_Preprocess_QC_config.yaml"

### Setup sample inputs
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"2_Reads/1_Untrimmed/*_1.fastq.gz")]

print("Detected the following samples:")
print(SAMPLE)

################################################################################
### Setup the desired outputs
rule all:
    input:
#        expand("3_Outputs/1_QC/2_CoverM/{sample}_coverM_mapped_host.tsv", sample=SAMPLE),
#        expand("3_Outputs/1_QC/3_nonpareil/{sample}.npo", sample=SAMPLE),
        "3_Outputs/1_QC/_report.tsv",
        "3_Outputs/1_QC/_nonpareil_metadata.tsv"
################################################################################
### Preprocess the reads using fastp
rule fastp:
    input:
        r1i = "2_Reads/1_Untrimmed/{sample}_1.fastq.gz",
        r2i = "2_Reads/1_Untrimmed/{sample}_2.fastq.gz"
    output:
        r1o = temp("2_Reads/2_Trimmed/{sample}_trimmed_1.fastq.gz"),
        r2o = temp("2_Reads/2_Trimmed/{sample}_trimmed_2.fastq.gz"),
        fastp_html = "2_Reads/3_fastp_results/{sample}.html",
        fastp_json = "2_Reads/3_fastp_results/{sample}.json"
    params:
        adapter1 = expand("{adapter1}", adapter1=config['adapter1']),
        adapter2 = expand("{adapter2}", adapter2=config['adapter2'])
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        8
    resources:
        mem_gb=24
    benchmark:
        "3_Outputs/0_Logs/{sample}_fastp.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_fastp.log"
    message:
        "Using FASTP to trim adapters and low quality sequences for {wildcards.sample}"
    shell:
        """
        fastp \
            --in1 {input.r1i} --in2 {input.r2i} \
            --out1 {output.r1o} --out2 {output.r2o} \
            --trim_poly_g \
            --trim_poly_x \
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
## Index host genomes:
rule index_ref:
    input:
        "1_References"
    output:
        bt2_index = "1_References/CattedRefs_renamed.fna.gz.rev.2.bt2l",
        rn_catted_ref = "1_References/CattedRefs_renamed.fna.gz"
    params:
        catted_ref = temp("1_References/CattedRefs.fna.gz")
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        40
    resources:
        mem_gb=180
    log:
        "3_Outputs/0_Logs/host_genome_indexing.log"
    message:
        "Concatenating and indexing host genomes with Bowtie2"
    shell:
        """
        # Concatenate input reference genomes
        cat {input}/*.gz > {params.catted_ref}

        # Add '_' separator for CoverM
        rename.sh in={params.catted_ref} \
        out={output.rn_catted_ref} \
        prefix=$(basename {input}) \
        -Xmx{resources.mem_gb}G

        # Index catted genomes
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {output.rn_catted_ref} {output.rn_catted_ref} \
        &> {log}
        """
################################################################################
### Map samples to host genomes, then split BAMs:
rule map_to_ref:
    input:
        r1i = "2_Reads/2_Trimmed/{sample}_trimmed_1.fastq.gz",
        r2i = "2_Reads/2_Trimmed/{sample}_trimmed_2.fastq.gz",
        catted_ref = "1_References/CattedRefs_renamed.fna.gz",
        bt2_index = "1_References/CattedRefs_renamed.fna.gz.rev.2.bt2l"
    output:
        all_bam = temp("3_Outputs/1_QC/1_BAMs/{sample}.bam"),
        host_bam = "3_Outputs/1_QC/1_Host_BAMs/{sample}_host.bam",
        non_host_r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq",
        non_host_r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq",
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        20
    resources:
        mem_gb=90
    benchmark:
        "3_Outputs/0_Logs/{sample}_mapping.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_mapping.log"
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
        non_host_r1 = "2_Reads/4_Host_removed/{sample}_M_1.fastq",
        non_host_r2 = "2_Reads/4_Host_removed/{sample}_M_2.fastq",
    output:
        npo = "3_Outputs/1_QC/3_nonpareil/{sample}.npo"
    params:
        sample = "3_Outputs/1_QC/3_nonpareil/{sample}",
        badsample_r1 = "2_Reads/5_Poor_samples/{sample}_M_1.fastq.gz",
        badsample_r2 = "2_Reads/5_Poor_samples/{sample}_M_2.fastq.gz"
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        10
    resources:
        mem_gb=45
    benchmark:
        "3_Outputs/0_Logs/{sample}_nonpareil.benchmark.tsv"
    message:
        "Estimating microbial diversity using nonpareil"
    shell:
        """
        mkdir -p 3_Outputs/1_QC/3_nonpareil

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
        if [ $(( $(stat -c '%s' {input.non_host_r1}) / 1024 / 1024 )) -lt 200 ]
        then
        mv {input.non_host_r1} {param.badsample_r1} && mv {input.non_host_r2} {param.badsample_r2}
        fi
        
        """
################################################################################
### Calculate % of each sample's reads mapping to host genome/s
rule coverM:
    input:
        bam = "3_Outputs/1_QC/1_BAMs/{sample}.bam",
        npo = "3_Outputs/1_QC/3_nonpareil/{sample}.npo"
    output:
        "3_Outputs/1_QC/2_CoverM/{sample}_coverM_mapped_host.tsv"
    params:
        assembly = "1_References/CattedRefs_renamed.fna.gz"
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        10
    resources:
        mem_gb=45
    benchmark:
        "3_Outputs/0_Logs/{sample}_coverM.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{sample}_coverM.log"
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
### Create summary table from outputs
rule report:
    input:
        coverm = expand("3_Outputs/1_QC/2_CoverM/{sample}_coverM_mapped_host.tsv", sample=SAMPLE),
        fastp = expand("2_Reads/3_fastp_results/{sample}.json", sample=SAMPLE)
    output:
        report = "3_Outputs/1_QC/_report.tsv",
        npar_metadata = "3_Outputs/1_QC/_nonpareil_metadata.tsv"
    params:
        tmpdir = "3_Outputs/1_QC/temp"
    conda:
        "1_Preprocess_QC.yaml"
    threads:
        1
    resources:
        mem_gb=45
    message:
        "Creating a final preprocessing report"
    shell:
        """
        #Create nonpareil sample metadata file
        mkdir -p {params.tmpdir}
        for i in {input.npar}; do echo $(basename $i) >> {params.tmpdir}/files.txt; done
        for i in {input.npar}; do echo $(basename ${{i/.npo/}}) >> {params.tmpdir}/names.txt; done
        for i in {input.npar}; do echo "#f03b20" >> {params.tmpdir}/colours.txt; done
        echo -e "File\tName\tColour" > {params.tmpdir}/headers.txt
        paste {params.tmpdir}/files.txt {params.tmpdir}/names.txt {params.tmpdir}/colours.txt > {params.tmpdir}/merged.tsv
        cat {params.tmpdir}/headers.txt {params.tmpdir}/merged.tsv > {output.npar_metadata}
        rm -r {params.tmpdir}

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

        paste {params.tmpdir}/names.tsv {params.tmpdir}/read_pre_filt.tsv {params.tmpdir}/read_post_filt.tsv {params.tmpdir}/bases_pre_filt.tsv {params.tmpdir}/bases_post_filt.tsv {params.tmpdir}/adapter_trimmed_reads.tsv {params.tmpdir}/adapter_trimmed_bases.tsv {params.tmpdir}/host_reads.tsv > {params.tmpdir}/preprocessing_stats.tsv
        echo -e "sample\treads_pre_filt\treads_post_filt\tbases_pre_filt\tbases_post_filt\tadapter_trimmed_reads\tadapter_trimmed_bases\thost_reads" > {params.tmpdir}/headers.tsv
        cat {params.tmpdir}/headers.tsv {params.tmpdir}/preprocessing_stats.tsv > {output.report}
        rm -r {params.tmpdir}
        """
################################################################################
onsuccess:
    shell("""
            mail -s "workflow completed" raph.eisenhofer@gmail.com < {log}
          """)
