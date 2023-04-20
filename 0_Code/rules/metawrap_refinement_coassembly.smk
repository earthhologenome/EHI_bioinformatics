################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        os.path.join(
            config["workdir"], 
            "{EHA}_binning/binning_complete"
            ),
    output:
        stats=os.path.join(
            config["workdir"],
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        contigmap=os.path.join(
            config["workdir"],
            "{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.contigs",
            )
    params:
        binning=os.path.join(config["workdir"] + "/{EHA}_binning"),
        outdir=os.path.join(config["workdir"] + "/{EHA}_refinement")
    threads: 16
    resources:
        mem_gb=128,
        time="06:00:00",
    benchmark:
        os.path.join(config["logdir"] + "/refinement_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/refinement_log_{EHA}.log")
    message:
        "Refining {wildcards.EHA} bins with MetaWRAP's bin refinement module"
    shell:
        """
        #Installing metawrap via conda is a pain in the arse, so using the module on Mjolnir here.
        module load metawrap-mg/1.3.2

        # Setup checkM path (needed for conda, not module)
        # export checkmdb={config[checkmdb]}
        # printf $checkmdb | checkm data setRoot
        
        metawrap bin_refinement \
            -m {resources.mem_gb} \
            -t {threads} \
            -o {params.outdir} \
            -A {params.binning}/concoct_bins/ \
            -B {params.binning}/maxbin2_bins/ \
            -C {params.binning}/metabat2_bins/ \
            -c 50 \
            -x 10

        # Rename metawrap bins to match coassembly group:
        mv {params.outdir}/metawrap_50_10_bins.stats {output.stats}
        mv {params.outdir}/metawrap_50_10_bins.contigs {output.contigmap}
        sed -i'' '2,$s/bin/{wildcards.EHA}_bin/g' {output.stats}
        sed -i'' 's/bin/{wildcards.EHA}_bin/g' {output.contigmap}
        for bin in {params.outdir}/metawrap_50_10_bins/*.fa;
            do mv $bin ${{bin/bin./{wildcards.EHA}_bin.}};
                done

        # Compress output bins
        pigz -p {threads} {params.outdir}/metawrap_50_10_bins/*.fa

        #Print the number of MAGs to a file for combining with the assembly report
        ls -l {params.outdir}/metawrap_50_10_bins/*.fa.gz | wc -l > {config[workdir]}/{wildcards.EHA}_bins.tsv

        # rm -r {params.binning}/work_files/
        # rm -r {params.outdir}/work_files/
        # rm {params.binning}/concoct_bins/*.fa
        # rm {params.binning}/maxbin2_bins/*.fa
        # rm {params.binning}/metabat2_bins/*.fa

        mv {params.outdir}/metawrap_50_10.stats {output.stats} 
        mv {params.outdir}/metawrap_50_10.contigs {output.contigmap}
        """