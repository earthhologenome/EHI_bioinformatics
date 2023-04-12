################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        os.path.join(
            config["workdir"], "{PRB}/", "{EHI}/", "{EHA}_binning/binning_complete"
        ),
    output:
        stats=os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_refinement/",
            "{EHA}_metawrap_70_10_bins.stats",
        ),
        contigmap=os.path.join(
            config["workdir"],
            "{PRB}/",
            "{EHI}/",
            "{EHA}_refinement/",
            "{EHA}_metawrap_70_10_bins.contigs",
        ),
    params:
        concoct="{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/concoct_bins",
        maxbin2="{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/maxbin2_bins",
        metabat2="{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/metabat2_bins",
        binning_wfs="{config['workdir']}/{PRB}/{EHI}/{EHA}_binning/work_files",
        refinement_wfs="{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/work_files",
        outdir="{config['workdir']}/{PRB}/{EHI}/{EHA}_refinement/",
    threads: 16
    resources:
        mem_gb=128,
        time="06:00:00",
    benchmark:
        "{{config['logdir']}}/binning_benchmark_{PRB}_{EHI}_{EHA}.tsv"
    log:
        "{{config['logdir']}}/binning_log_{PRB}_{EHI}_{EHA}.log",
    message:
        "Refining {wildcards.EHA} bins with MetaWRAP's bin refinement module"
    shell:
        """
        #Installing metawrap via conda is a pain in the arse, so using the module on Mjolnir here.
        module load metawrap-mg/1.3.2


        # Setup checkM path
        export checkmdb={config[checkmdb]}
        printf $checkmdb | checkm data setRoot
        
        metawrap bin_refinement \
            -m {resources.mem_gb} \
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
        sed -i'' '2,$s/bin/{wildcards.EHA}_bin/g' {output.stats}
        sed -i'' 's/bin/{wildcards.EHA}_bin/g' {output.contigmap}
        for bin in {params.outdir}/metawrap_70_10_bins/*.fa;
            do mv $bin ${{bin/bin./{wildcards.EHA}_bin.}};
                done

        # Compress output bins
        pigz -p {threads} {params.outdir}/*bins/*.fa

        rm -r {params.binning_wfs}
        rm -r {params.refinement_wfs}
        rm {params.concoct}/*.fa
        rm {params.maxbin2}/*.fa
        rm {params.metabat2}/*.fa
        """