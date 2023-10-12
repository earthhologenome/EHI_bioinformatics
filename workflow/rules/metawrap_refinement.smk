################################################################################
### Automatically refine bins using metaWRAP's refinement module
rule metaWRAP_refinement:
    input:
        concoct=os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_{EHA}_binning/concoct_binning_complete"
            ),
        metabat2=os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_{EHA}_binning/metabat2_binning_complete"
            ),
        maxbin2=os.path.join(
            config["workdir"], 
            "{PRB}_{EHI}_{EHA}_binning/maxbin2_binning_complete"
            ),
        contigs=os.path.join(
            config["workdir"], "{PRB}_{EHI}_assembly/", "{EHA}_contigs.fasta"
            )
    output:
        stats=os.path.join(
            config["workdir"],
            "{PRB}_{EHI}_{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.stats",
            ),
        contigmap=os.path.join(
            config["workdir"],
            "{PRB}_{EHI}_{EHA}_refinement/",
            "{EHA}_metawrap_50_10_bins.contigs",
            )
    params:
        binning=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_binning"),
        outdir=os.path.join(config["workdir"] + "/{PRB}_{EHI}_{EHA}_refinement"),
        stats_dir=directory(os.path.join(config["workdir"], "{EHA}_stats/")),
        contigsize=config["contigsize"]
    conda:
        f"{config['codedir']}/conda_envs/checkm2.yaml"
    threads: 16
    resources:
        mem_gb=80,
        time=estimate_time_refinement,
    benchmark:
        os.path.join(config["logdir"] + "/refinement_benchmark_{PRB}_{EHI}_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/refinement_log_{PRB}_{EHI}_{EHA}.log")
    message:
        "Refining {wildcards.EHA} bins with MetaWRAP's bin refinement module"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt {params.contigsize} ]
        then
            touch {output}
            
        elif
            [ $(ls -l {params.binning}/metabat2_bins/*.fa | wc -l) -lt 6 ]

        then
            touch {output}
        
        else
        
            # setup checkm2 db path (in case of first run)
            checkm2 database --setdblocation {config[checkmdb]}

            bash {config[codedir]}/scripts/metawrap_refinement.sh \
                -t {threads} \
                -o {params.outdir} \
                -A {params.binning}/concoct_bins/ \
                -B {params.binning}/maxbin2_bins/ \
                -C {params.binning}/metabat2_bins/ \
                -c 50 \
                -x 10

            # Rename output files, and sort metawrap by bin name
            head -1 {params.outdir}/metawrap_50_10_bins.stats > {params.outdir}/mw_colnames.tsv
            sed '1d;' {params.outdir}/metawrap_50_10_bins.stats | sort -k1,1 -t$'\t' > {params.outdir}/mw_sorted.tsv
            cat {params.outdir}/mw_colnames.tsv {params.outdir}/mw_sorted.tsv > {params.outdir}/mw_sorted_col.tsv
            mv {params.outdir}/mw_sorted_col.tsv {output.stats}
            mv {params.outdir}/metawrap_50_10_bins.contigs {output.contigmap}
            sed -i'' '2,$s/bin/{wildcards.EHA}_bin/g' {output.stats}
            sed -i'' 's/bin/{wildcards.EHA}_bin/g' {output.contigmap}

            # Rename metawrap bins to match coassembly group:
            for bin in {params.outdir}/metawrap_50_10_bins/*.fa;
                do mv $bin ${{bin/bin./{wildcards.EHA}_bin.}};
            done

            # Compress output bins
            pigz -p {threads} {params.outdir}/metawrap_50_10_bins/*.fa

            #Print the number of MAGs to a file for combining with the assembly report
            mkdir -p {params.stats_dir}
            ls -l {params.outdir}/metawrap_50_10_bins/*.fa.gz | wc -l > {params.stats_dir}/{wildcards.EHA}_bins.tsv

            # Reformat MAG headers for CoverM
            for mag in {params.outdir}/metawrap_50_10_bins/*.fa.gz;
                do rename.sh \
                    in=$mag \
                    out={params.outdir}/$(basename ${{mag/.fa.gz/_renamed.fa.gz}}) \
                    zl=9 \
                    prefix=$(basename ${{mag/.fa.gz/^}});
            done

            rm {params.outdir}/metawrap_50_10_bins/*.fa.gz
            for mag in {params.outdir}/*.fa.gz;
                do mv $mag {params.outdir}/metawrap_50_10_bins/$(basename ${{mag/_renamed/}});
            done

            # rm -r {params.binning}/work_files/
            # rm -r {params.outdir}/work_files/
            # rm {params.binning}/concoct_bins/*.fa
            # rm {params.binning}/maxbin2_bins/*.fa
            # rm {params.binning}/metabat2_bins/*.fa
        fi
        """