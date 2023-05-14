rule check_assembly_size:
    input:
        os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta")
    output:
        os.path.join(config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_checkpoint")
    threads: 1
    resources:
        mem_gb=8,
        time="00:05:00",
    shell:
        """
        # check the size of the output file
        if [ $(( $(stat -c '%s' {input}) / 1024 / 1024 )) -lt 50 ]
        then
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_uploaded
            mkdir -p {config[workdir]}/bams
            touch {config[workdir]}/bams/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}.bam
            mkdir -p {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_binning
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_binning/binning_complete
            mkdir -p {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/
            touch {config[workdir]}/{wildcards.PRB}/{wildcards.EHI}/{wildcards.EHA}_gtdbtk_combined_summary.tsv
            mkdir -p {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_refinement
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_refinement/{wildcards.EHA}_metawrap_50_10_bins.stats            
            mkdir -p {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_assembly/
            touch {config[workdir]}/{wildcards.PRB}_{wildcards.EHI}_assembly/{wildcards.EHA}_contigs.fasta.gz
     


            # Create sample table for airtable
            echo -e "sample\tEHA_number\tEHI_number\tN50\tL50\tnum_contigs\tlargest_contig\tassembly_length\tnum_bins\tassembly_mapping_percent" > headers.tsv

            #Add sample IDs
            echo {wildcards.EHI}_{wildcards.EHA} >> {wildcards.EHA}_sample_ids.tsv
            echo {wildcards.EHA} >> {wildcards.EHA}_EHA_ids.tsv
            echo {wildcards.EHI} >> {wildcards.EHA}_EHI_ids.tsv

            echo -e "0\t0\t0\t0\t0\t0\t0" >> {wildcards.EHA}_null.tsv

            paste {wildcards.EHA}_sample_ids.tsv {wildcards.EHA}_EHA_ids.tsv {wildcards.EHA}_EHI_ids.tsv {wildcards.EHA}_null.tsv > {wildcards.EHA}_nullz.tsv
            cat headers.tsv {wildcards.EHA}_nullz.tsv > /projects/ehi/data/REP/{wildcards.PRB}_{wildcards.EHI}_{wildcards.EHA}_final_stats.tsv

            touch {output}
            echo "Output file size is below minimum threshold. Stopping processing of this sample." >&2
            exit 1
        else
            touch {output}
        fi
        """
onerrmsg = "Processing of sample {wildcards.EHI} stopped because the assembly size is below the minimum threshold"
onerror:
    raise RuntimeError(onerrmsg)