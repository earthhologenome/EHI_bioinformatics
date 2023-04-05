################################################################################
################################################################################
################################################################################
# EHI snakefile for assembly/binning and MAG annotation
# Raphael Eisenhofer 04/2023
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


### Setup sample inputs, config, and working directory

configfile: "2_Assembly_Binning_config.yaml"

import pandas as pd

## The input will be automatically generated prior to the snakefile being launched-
## using the 'XXXXX.py' script, which pulls the information from AirTable and saves-
## it as 'abb_input.tsv'.

samples = pd.read_csv('abb_input.tsv', sep='\t')
valid_combinations = [(row['ID'], row['EHI_number']) for _, row in samples.iterrows()]

################################################################################
### Setup the desired outputs
rule all:
    params:
        eha_samples = valid_combinations
    input:
        expand(f"{config['workdir']}/{eha}/{sample}_1.fq.gz", 
               sample=valid_combinations[i][1], eha=valid_combinations[i][0], 
               i=range(len(valid_combinations)))
################################################################################
### Create EHA folder on ERDA
rule create_ASB_folder:
    output:
        f"{config['workdir']}/{config['abb']}/ERDA_folder_created"
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time='00:05:00'
    message:
        "Creating assembly batch folder on ERDA"
    shell:
        """
        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/ASB/{config['abb']} ; bye"
        touch {output}

        #Also, log the AirTable that the ASB is running!
        python {config['codedir']}/airtable/log_asb_start_airtable.py --code={config['abb']}

        """
################################################################################
### Fetch preprocessed reads from ERDA
rule download_from_ERDA:
    input:
        f"{config['workdir']}/{config['abb']}/ERDA_folder_created"
    output:
        r1 = f"{config['workdir']}/{eha}/{sample}_1.fq.gz",
        r2 = f"{config['workdir']}/{eha}/{sample}_2.fq.gz",
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time='00:15:00'
    message:
        "Fetching {wildcards.sample} from ERDA"
    shell:
        """
        mkdir -p {params.eha}
        
        lftp sftp://erda -e "mirror --include-glob='{wildcards.prbatch}/{sample}*.fq.gz' /EarthHologenomeInitiative/Data/RAW/ {params.eha_folder}; bye"

        """



# ################################################################################
# ### Perform Coassemblies on each sample group
# rule Coassembly:
#     input:
#         reads = "2_Reads/4_Host_removed/{group}"
#     output:
#         Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
#     params:
#         workdir = "3_Outputs/2_Coassemblies/{group}",
#         r1_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_1.fq.gz"),
#         r2_cat = temp("3_Outputs/2_Coassemblies/{group}/{group}_2.fq.gz"),
#         assembler = expand("{assembler}", assembler=config['assembler']),
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     threads:
#         24
#     resources:
#         mem_gb=196,
#         time='36:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_coassembly.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_coassembly.log"
#     message:
#         "Coassembling {wildcards.group} using {params.assembler}"
#     shell:
#         """
#         # Set up input reads variable for megahit
#         R1=$(for i in 2_Reads/4_Host_removed/{wildcards.group}/*_1.fq.gz; do echo $i | tr '\n' ,; done)
#         R2=$(for i in 2_Reads/4_Host_removed/{wildcards.group}/*_2.fq.gz; do echo $i | tr '\n' ,; done)

#         # Run megahit
#             megahit \
#                 -t {threads} \
#                 --verbose \
#                 --min-contig-len 1500 \
#                 -1 $R1 -2 $R2 \
#                 -f \
#                 -o {params.workdir}
#                 2> {log}

#         # Move the Coassembly to final destination
#             mv {params.workdir}/final.contigs.fa {output.Coassembly}

#         # Reformat headers
#             sed -i 's/ /-/g' {output.Coassembly}

#         """
# ################################################################################
# ### Create QUAST reports of coassemblies
# rule QUAST:
#     input:
#         Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
#     output:
#         report = directory("3_Outputs/2_Coassemblies/{group}_QUAST"),
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     threads:
#         8
#     resources:
#         mem_gb=48,
#         time='01:00:00'
#     message:
#         "Running -QUAST on {wildcards.group} coassembly"
#     shell:
#         """
#         # Run QUAST
#         quast \
#             -o {output.report} \
#             --threads {threads} \
#             {input.Coassembly}

#         # # Rename QUAST files
#         # for i in {output.report}/*;
#         #     do mv $i {output.report}/{wildcards.group}_$(basename $i);
#         #         done
        
#         # Parse select metrics for final report
#         grep N50 {output.report}/report.tsv | cut -f2 > {output.report}/n50.tsv
#         grep L50 {output.report}/report.tsv | cut -f2 > {output.report}/l50.tsv
#         grep "# contigs (>= 0 bp)" {output.report}/report.tsv | cut -f2 > {output.report}/ncontigs.tsv
#         grep "Largest contig" {output.report}/report.tsv | cut -f2 > {output.report}/largestcontig.tsv
#         grep "Total length (>= 0 bp)" {output.report}/report.tsv | cut -f2 > {output.report}/totallength.tsv

#         # paste into a single table
#         paste {output.report}/n50.tsv \
#               {output.report}/l50.tsv \
#               {output.report}/ncontigs.tsv \
#               {output.report}/largestcontig.tsv \
#               {output.report}/totallength.tsv > {output.report}/{wildcards.group}_assembly_report.tsv
#         """
# ################################################################################
# ### Map reads to the coassemblies
# rule Coassembly_index:
#     input:
#         report = "3_Outputs/2_Coassemblies/{group}_QUAST/"
#     output:
#         bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l",
#     params:
#         Coassembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta"
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     threads:
#         24
#     resources:
#         mem_gb=128,
#         time='04:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_coassembly_indexing.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_coassembly_indexing.log"
#     message:
#         "Indexing {wildcards.group} coassembly using Bowtie2"
#     shell:
#         """
#         # Index the coassembly
#         bowtie2-build \
#             --large-index \
#             --threads {threads} \
#             {params.Coassembly} {params.Coassembly} \
#         &> {log}
#         """
# ################################################################################
# ### Map reads to the coassemblies
# rule Coassembly_mapping:
#     input:
#         bt2_index = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta.rev.2.bt2l"
#     output:
#         bam = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}/{sample}.bam"  
#     params:
#         r1 = "2_Reads/4_Host_removed/{group}/{sample}_M_1.fq.gz",
#         r2 = "2_Reads/4_Host_removed/{group}/{sample}_M_2.fq.gz",
#         assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     threads:
#         8
#     resources:
#         mem_gb=40,
#         time='08:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_{sample}_coassembly_mapping.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_{sample}_coassembly_mapping.log"
#     message:
#         "Mapping {wildcards.sample} to {wildcards.group} coassembly using Bowtie2"
#     shell:
#         """
#         # Map reads to catted reference using Bowtie2
#         bowtie2 \
#             --time \
#             --threads {threads} \
#             -x {params.assembly} \
#             -1 {params.r1} \
#             -2 {params.r2} \
#         | samtools sort -@ {threads} -o {output.bam}
#         """
# ################################################################################
# ### Bin contigs using metaWRAP's binning module
# rule metaWRAP_binning:
#     input:
#         bams = lambda wildcards: ["3_Outputs/3_Coassembly_Mapping/BAMs/{}/{}.bam/".format(wildcards.group, sample) for sample in GROUPS[wildcards.group]]
#     output:
#         "3_Outputs/4_Binning/{group}/Done.txt"
#     params:
#         concoct = "3_Outputs/4_Binning/{group}/concoct_bins",
#         maxbin2 = "3_Outputs/4_Binning/{group}/maxbin2_bins",
#         metabat2 = "3_Outputs/4_Binning/{group}/metabat2_bins",
#         outdir = "3_Outputs/4_Binning/{group}",
#         bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
#         assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
#         memory = "180"
#     conda:
#         "conda_envs/2_MetaWRAP.yaml"
#     threads:
#         32
#     resources:
#         mem_gb=256,
#         time='96:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_coassembly_binning.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_coassembly_binning.log"
#     message:
#         "Binning {wildcards.group} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
#     shell:
#         """
#         # Create dummy fq/assembly files to trick metaWRAP into running without mapping
#         mkdir -p {params.outdir}/work_files

#         touch {params.outdir}/work_files/assembly.fa.bwt

#         for bam in {params.bams}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
#         for bam in {params.bams}/*.bam; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done


#         #Symlink BAMs for metaWRAP
#         for bam in {params.bams}/*.bam; do ln -sf `pwd`/$bam {params.outdir}/work_files/$(basename $bam); done

#         # Run metaWRAP binning
#         metawrap binning -o {params.outdir} \
#             -t {threads} \
#             -m {params.memory} \
#             -a {params.assembly} \
#             -l 1500 \
#             --metabat2 \
#             --maxbin2 \
#             --concoct \
#         {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq

#         # Create dummy file for refinement input
#         echo "Binning complete" > {output}
#         """
# ################################################################################
# ### Automatically refine bins using metaWRAP's refinement module
# rule metaWRAP_refinement:
#     input:
#         "3_Outputs/4_Binning/{group}/Done.txt"
#     output:
#         stats = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats",
#         contigmap = "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.contigs"
#     params:
#         concoct = "3_Outputs/4_Binning/{group}/concoct_bins",
#         maxbin2 = "3_Outputs/4_Binning/{group}/maxbin2_bins",
#         metabat2 = "3_Outputs/4_Binning/{group}/metabat2_bins",
#         binning_wfs = "3_Outputs/4_Binning/{group}/work_files",
#         refinement_wfs = "3_Outputs/5_Refined_Bins/{group}/work_files",
#         outdir = "3_Outputs/5_Refined_Bins/{group}",
#         memory = "256",
#         group = "{group}"
#     conda:
#         "conda_envs/2_MetaWRAP.yaml"
#     threads:
#         32
#     resources:
#         mem_gb=256,
#         time='36:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
#     message:
#         "Refining {wildcards.group} bins with MetaWRAP's bin refinement module"
#     shell:
#         """
#         # Setup checkM path
#         export checkmdb={config[checkmdb]}
#         printf $checkmdb | checkm data setRoot
        
#         metawrap bin_refinement \
#             -m {params.memory} \
#             -t {threads} \
#             -o {params.outdir} \
#             -A {params.concoct} \
#             -B {params.maxbin2} \
#             -C {params.metabat2} \
#             -c 70 \
#             -x 10

#         # Rename metawrap bins to match coassembly group:
#         mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
#         mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
#         sed -i'' '2,$s/bin/{params.group}_bin/g' {output.stats}
#         sed -i'' 's/bin/{params.group}_bin/g' {output.contigmap}
#         for bin in {params.outdir}/metawrap_70_10_bins/*.fa;
#             do mv $bin ${{bin/bin./{params.group}_bin.}};
#                 done

#         # Compress output bins
#         pigz -p {threads} {params.outdir}/*bins/*.fa

#         rm -r {params.binning_wfs}
#         rm -r {params.refinement_wfs}
#         rm {params.concoct}/*.fa
#         rm {params.maxbin2}/*.fa
#         rm {params.metabat2}/*.fa
#         """
# ################################################################################
# ### Calculate the number of reads that mapped to coassemblies
# rule coverM_assembly:
#     input:
#         "3_Outputs/5_Refined_Bins/{group}/{group}_metawrap_70_10_bins.stats"
#     output:
#         coverm = "3_Outputs/6_Coassembly_CoverM/{group}_assembly_coverM.txt",
#         euk = "3_Outputs/6_Coassembly_CoverM/{group}_eukaryotic_coverM.tsv"
#     params:
#         mapped_bams = "3_Outputs/3_Coassembly_Mapping/BAMs/{group}",
#         assembly = "3_Outputs/2_Coassemblies/{group}/{group}_contigs.fasta",
#         binning_files = "3_Outputs/4_Binning/{group}",
#         refinement_files = "3_Outputs/5_Refined_Bins/{group}",
#         group = "{group}"
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     threads:
#         24
#     resources:
#         mem_gb=128,
#         time='06:00:00'
#     benchmark:
#         "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
#     log:
#         "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
#     message:
#         "Calculating coassembly mapping rate for {wildcards.group} with CoverM"
#     shell:
#         """
#         coverm genome \
#             -b {params.mapped_bams}/*.bam \
#             --genome-fasta-files {params.assembly} \
#             -m relative_abundance \
#             -t {threads} \
#             --min-covered-fraction 0 \
#             > {output.coverm}

#         #Run coverm for the eukaryotic assessment pipeline
#         coverm genome \
#             -s - \
#             -b {params.mapped_bams}/*.bam \
#             -m relative_abundance count mean covered_fraction \
#             -t {threads} \
#             --min-covered-fraction 0 \
#             > {output.euk}

#         # Create directory for dereplication groups:
#         mkdir -p 3_Outputs/5_Refined_Bins/dRep_groups

#         #Print the number of MAGs to a file for combining with the assembly report
#         ls -l {params.refinement_files}/metawrap_70_10_bins/*.fa.gz | wc -l > {params.group}_bins.tsv;
#         """
# ################################################################################
# ### Generate output summary table
# rule generate_summary:
#     input:
#         "3_Outputs/6_Coassembly_CoverM/{group}_assembly_coverM.txt"
#     output:
#         "3_Outputs/{group}_coassembly_summary.tsv"
#     conda:
#         "conda_envs/2_Assembly_Binning.yaml"
#     params:
#         group = "{group}"
#     threads:
#         1
#     resources:
#         mem_gb=16,
#         time='00:05:00'
#     message:
#         "Creating final coassembly summary table"
#     shell:
#         """
#         #Create the final output summary table
#         #parse QUAST outputs for assembly stats
#         echo -e "sample\tN50\tL50\tnum_contigs\tlargest_contig\ttotal_length\tnum_bins\tassembly_mapping_percent" > headers.tsv

     
#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do cat 3_Outputs/2_Coassemblies/{params.group}_QUAST/{params.group}_assembly_report.tsv >> {params.group}_temp_report.tsv;
#         done


#         #Add in the % mapping to assembly stats
#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do echo $(basename ${{sample/.bam/}}) >> {params.group}_sample_ids.tsv;
#         done

#         paste {params.group}_sample_ids.tsv {params.group}_temp_report.tsv > {params.group}_temp2_report.tsv

#         #Add in the # of bins
#         for sample in 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam;
#             do cat {params.group}_bins.tsv >> {params.group}_number_bins.tsv;
#         done

#         paste {params.group}_temp2_report.tsv {params.group}_number_bins.tsv > {params.group}_temp3_report.tsv

#         ls -l 3_Outputs/3_Coassembly_Mapping/BAMs/{params.group}/*.bam | wc -l > {params.group}_n_samples.tsv

#         nsamples=$( cat {params.group}_n_samples.tsv )
#         nsamples1=$(( nsamples + 1 ))
#         echo $nsamples1
#         for sample in `seq 2 $nsamples1`;
#             do cut -f"$sample" 3_Outputs/6_Coassembly_CoverM/{params.group}_assembly_coverM.txt | sed -n 3p >> {params.group}_relabun.tsv;
#         done

#         paste {params.group}_temp3_report.tsv {params.group}_relabun.tsv > {params.group}_temp4_report.tsv

#         #Combine them into the final assembly report
#         cat headers.tsv {params.group}_temp4_report.tsv > {output}
#         """
# ################################################################################
# ### Clean up
# rule clean:
#     input:
#         expand("3_Outputs/{group}_coassembly_summary.tsv", group=GROUPS.keys())
#     output:
#         "3_Outputs/pipeline_complete.txt"
#     threads:
#         1
#     resources:
#         mem_gb=16,
#         time='00:05:00'
#     message:
#         "Cleaning up temp files"
#     shell:
#         """
#         rm *.tsv
#         touch {output}
#         """