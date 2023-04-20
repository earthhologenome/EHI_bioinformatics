################################################################################
### Create QUAST reports of coassemblies
rule QUAST:
    input:
        os.path.join(config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"),
    output:
        directory(os.path.join(config["workdir"], "{EHA}_QUAST")),
    conda:
        f"{config['codedir']}/conda_envs/assembly_binning.yaml"
    threads: 4
    resources:
        mem_gb=32,
        time="00:30:00",
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
              {output}/totallength.tsv > {output}/{wildcards.EHA}_assembly_report.tsv
        """