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

df = pd.read_csv('abb_input.tsv', sep='\t')

# Use set to create a list of valid combinations of wildcards. Note that 'ID' = EHA number.
valid_combinations = set((row['PR_batch'], row['EHI_number'], row['ID']) for _, row in df.iterrows())


################################################################################
### Setup the desired outputs
rule all:
    input:
        expand("downloaded_data/{combo[0]}/{combo[1]}/{combo[2]}_contigs.fasta.gz", combo=valid_combinations)

################################################################################
### Create EHA folder on ERDA
rule create_ASB_folder:
    output:
        f"{config['workdir']}/ERDA_folder_created"
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
        lftp sftp://erda -e "mkdir -f EarthHologenomeInitiative/Data/ASB/{{config['abb']}} ; bye"
        touch {output}

        #Also, log the AirTable that the ASB is running!
        python {{config['codedir']}}/airtable/log_asb_start_airtable.py --code={{config['abb']}}
        """
################################################################################
### Fetch preprocessed reads from ERDA
rule download_from_ERDA:
    input:
        f"{config['workdir']}/ERDA_folder_created"
    output:
        r1 = f"{config['workdir']}/{PRB}/{EHI}_M_1.fq.gz",
        r2 = f"{config['workdir']}/{PRB}/{EHI}_M_2.fq.gz"
    conda:
        f"{config['codedir']}/conda_envs/lftp.yaml"
    threads:
        1
    resources:
        load=8,
        mem_gb=8,
        time='00:15:00'
    message:
        "Fetching metagenomics reads for {wildcards.EHI} from ERDA"
    shell:
        """
        lftp sftp://erda -e "mirror --include-glob='{wildcards.PRB}/{wildcards.EHI}*.fq.gz' /EarthHologenomeInitiative/Data/PPR/ {{config['workdir']}}/; bye"
        """

################################################################################
### Perform individual assembly on each sample
rule assembly:
    input:
        r1 = f"{config['workdir']}/{PRB}/{EHI}_M_1.fq.gz",
        r2 = f"{config['workdir']}/{PRB}/{EHI}_M_2.fq.gz"
    output:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta"
    params:
        assembler = expand("{assembler}", assembler=config['assembler']),
    conda:
        f"{config['codedir']}/conda_envs/2_Assembly_Binning_config.yaml"
    threads:
        16
    resources:
        mem_gb=128,
        time='36:00:00'
    benchmark:
        f"{config['logdir']}/assembly_benchmark_{EHA}.tsv"
    log:
        f"{config['logdir']}/assembly_log_{EHA}.log"
    message:
        "Assembling {wildcards.EHI} using {params.assembler}"
    shell:
        """
        # Run megahit
            megahit \
                -t {threads} \
                --verbose \
                --min-contig-len 1500 \
                -1 {input.r1} -2 {input.r2} \
                -f \
                -o {{config['workdir']}}/{PRB}/{EHI}
                2> {log}

        # Move the Coassembly to final destination
            mv {{config['workdir']}}/{PRB}/{EHI}/final.contigs.fa {output}

        # Reformat headers
            sed -i 's/ /-/g' {output}
        """
################################################################################
### Create QUAST reports of coassemblies
rule QUAST:
    input:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta"
    output:
        directory(f"{config['workdir']}/{PRB}/{EHI}/{EHA}_QUAST"),
    conda:
        f"{config['codedir']}/conda_envs/2_Assembly_Binning_config.yaml"
    threads:
        4
    resources:
        mem_gb=32,
        time='00:30:00'
    message:
        "Running -QUAST on {wildcards.EHA} coassembly"
    shell:
        """
        # Run QUAST
        quast \
            -o {output} \
            --threads {threads} \
            {input}
      
        # Parse select metrics for final report
        grep N50 {output}/report.tsv | cut -f2 > {output}/n50.tsv
        grep L50 {output}/report.tsv | cut -f2 > {output}/l50.tsv
        grep "# contigs (>= 0 bp)" {output}/report.tsv | cut -f2 > {output}/ncontigs.tsv
        grep "Largest contig" {output}/report.tsv | cut -f2 > {output}/largestcontig.tsv
        grep "Total length (>= 0 bp)" {output}/report.tsv | cut -f2 > {output}/totallength.tsv

        # paste into a single table
        paste {output}/n50.tsv \
              {output}/l50.tsv \
              {output}/ncontigs.tsv \
              {output}/largestcontig.tsv \
              {output}/totallength.tsv > {output}/{wildcards.group}_assembly_report.tsv
        """
################################################################################
### Index assemblies
rule assembly_index:
    input:
        directory(f"{config['workdir']}/{PRB}/{EHI}/{EHA}_QUAST")
    output:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta.rev.2.bt2l"
    params:
        contigs = f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta"
    conda:
        f"{config['codedir']}/conda_envs/2_Assembly_Binning_config.yaml"
    threads:
        16
    resources:
        mem_gb=96,
        time='02:00:00'
    benchmark:
        f"{config['logdir']}/index_assembly_benchmark_{EHA}.tsv"
    log:
        f"{config['logdir']}/index_assembly_log_{EHA}.log"
    message:
        "Indexing {wildcards.EHA} assembly using Bowtie2"
    shell:
        """
        # Index the coassembly
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {params.contigs} {params.contigs} \
        &> {log}
        """
################################################################################
### Map reads to assemblies
rule assembly_mapping:
    input:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta.rev.2.bt2l",
        r1 = f"{config['workdir']}/{PRB}/{EHI}_M_1.fq.gz",
        r2 = f"{config['workdir']}/{PRB}/{EHI}_M_2.fq.gz",
        contigs = f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta"
    output:
        f"{config['workdir']}/{PRB}/{EHI}/{EHI}_{EHA}.bam"
    conda:
        f"{config['codedir']}/conda_envs/2_Assembly_Binning_config.yaml"
    threads:
        16
    resources:
        mem_gb=48,
        time='05:00:00'
    benchmark:
        f"{config['logdir']}/assembly_mapping_benchmark_{EHA}.tsv"
    log:
        f"{config['logdir']}/assembly_mapping_log_{EHA}.log"
    message:
        "Mapping {wildcards.EHI} to {wildcards.EHA} assembly using Bowtie2"
    shell:
        """
        # Map reads to assembly using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.assembly} \
            -1 {input.r1} \
            -2 {input.r2} \
        | samtools sort -@ {threads} -o {output}
        """
################################################################################
### Bin contigs using metaWRAP's binning module
rule metaWRAP_binning:
    input:
        bam = f"{config['workdir']}/{PRB}/{EHI}/{EHI}_{EHA}.bam",
        contigs = f"{config['workdir']}/{PRB}/{EHI}/{EHA}_contigs.fasta"
    output:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/binning_complete"
    params:
        outdir = directory(f"{config['workdir']}/{PRB}/{EHI}/{EHA}_binning")
    conda:
        f"{config['codedir']}/conda_envs/2_MetaWRAP.yaml"
    threads:
        16
    resources:
        mem_gb=96,
        time='06:00:00'
    benchmark:
        f"{config['logdir']}/binning_benchmark_{EHA}.tsv"
    log:
        f"{config['logdir']}/binning_log_{EHA}.log"
    message:
        "Binning {wildcards.EHA} contigs with MetaWRAP (concoct, maxbin2, metabat2)"
    shell:
        """
        # Create dummy fq/assembly files to trick metaWRAP into running without mapping
        mkdir -p {params.outdir}
        mkdir -p {params.outdir}/work_files

        touch {params.outdir}/work_files/assembly.fa.bwt

        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_1.fastq}}); done
        for bam in {input.bam}; do echo "@" > {params.outdir}/work_files/$(basename ${{bam/.bam/_2.fastq}}); done

        #Symlink BAMs for metaWRAP
        for bam in {input.bam}; do ln -sf $bam {params.outdir}/work_files/$(basename $bam); done

        # Run metaWRAP binning
        metawrap binning -o {params.outdir} \
            -t {threads} \
            -m {resources.mem_gb} \
            -a {input.contigs} \
            -l 1500 \
            --metabat2 \
            --maxbin2 \
            --concoct \
        {params.outdir}/work_files/*_1.fastq {params.outdir}/work_files/*_2.fastq

        # Create output for the next rule
        touch {output}
        """
################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        f"{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/binning_complete"
    output:
        stats = f"{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/{EHA}_metawrap_70_10_bins.stats",
        contigmap = f"{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/{EHA}_metawrap_70_10_bins.contigs"
    params:
        concoct = "3_Outputs/4_Binning/{group}/concoct_bins",
        maxbin2 = "3_Outputs/4_Binning/{group}/maxbin2_bins",
        metabat2 = "3_Outputs/4_Binning/{group}/metabat2_bins",
        binning_wfs = "3_Outputs/4_Binning/{group}/work_files",
        refinement_wfs = "3_Outputs/5_Refined_Bins/{group}/work_files",
        outdir = "3_Outputs/5_Refined_Bins/{group}",
        memory = "256",
        group = "{group}"
    conda:
        "conda_envs/2_MetaWRAP.yaml"
    threads:
        32
    resources:
        mem_gb=256,
        time='36:00:00'
    benchmark:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.benchmark.tsv"
    log:
        "3_Outputs/0_Logs/{group}_coassembly_bin_refinement.log"
    message:
        "Refining {wildcards.group} bins with MetaWRAP's bin refinement module"
    shell:
        """
        # Setup checkM path
        export checkmdb={config[checkmdb]}
        printf $checkmdb | checkm data setRoot
        
        metawrap bin_refinement \
            -m {params.memory} \
            -t {threads} \
            -o {params.outdir} \
            -A {params.concoct} \
            -B {params.maxbin2} \
            -C {params.metabat2} \
            -c 70 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        mv {params.outdir}/metawrap_70_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_70_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/{params.group}_bin/g' {output.stats}
        sed -i'' 's/bin/{params.group}_bin/g' {output.contigmap}
        for bin in {params.outdir}/metawrap_70_10_bins/*.fa;
            do mv $bin ${{bin/bin./{params.group}_bin.}};
                done

        # Compress output bins
        pigz -p {threads} {params.outdir}/*bins/*.fa

        rm -r {params.binning_wfs}
        rm -r {params.refinement_wfs}
        rm {params.concoct}/*.fa
        rm {params.maxbin2}/*.fa
        rm {params.metabat2}/*.fa
        """
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