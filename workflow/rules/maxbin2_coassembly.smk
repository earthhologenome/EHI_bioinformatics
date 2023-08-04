################################################################################
### Bin contigs using maxbin2
### Note that I've copied how the MetaWRAP pipeline calls MaxBin2 (https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-modules/binning.sh)
rule maxbin2:
    input:
        bam=expand(os.path.join(
            config["workdir"], "bams/", "{combo[0]}_{combo[1]}_{combo[2]}.bam"
            ),
            combo=valid_combinations
        ),
        contigs=os.path.join(
            config["workdir"], "{EHA}_assembly/", "{EHA}_contigs.fasta"
            ),
    output:
        os.path.join(
            config["workdir"], "{EHA}_binning/maxbin2_binning_complete"
            ),
    params:
        outdir=os.path.join(config["workdir"] + "/{EHA}_binning"),
        contigsize=config["contigsize"]
    conda:
        f"{config['codedir']}/conda_envs/maxbin2.yaml"
    threads: 16
    resources:
        mem_gb=64,
        time=estimate_time_binning,
    benchmark:
        os.path.join(config["logdir"] + "/binning_benchmark_{EHA}.tsv")
    log:
        os.path.join(config["logdir"] + "/binning_log_{EHA}.log")
    message:
        "Binning {wildcards.EHA} contigs with maxbin2"
    shell:
        """
        if [ $(( $(stat -c '%s' {input.contigs}) / 1024 / 1024 )) -lt {params.contigsize} ]
        then
            touch {output}

        else

            # summarise contig depths
            jgi_summarize_bam_contig_depths \
            --outputDepth {params.outdir}/mb2_master_depth.txt \
            --noIntraDepthVariance \
            {input.bam}

            #calculate total numper of columns
            A=($(head -n 1 {params.outdir}/mb2_master_depth.txt)) 
            N=${#A[*]}

            # split the contig depth file into multiple files
            if [ -f {params.outdir}/mb2_abund_list.txt ]; then rm {params.outdir}/mb2_abund_list.txt; fi
            for i in $(seq 4 $N); do 
                sample=$(head -n 1 {params.outdir}/mb2_master_depth.txt | cut -f $i)
                grep -v totalAvgDepth {params.outdir}/mb2_master_depth.txt | cut -f 1,$i > {params.outdir}/mb2_${sample%.*}.txt
                if [[ $out == /* ]]; then
                    echo {params.outdir}/mb2_${sample%.*}.txt >> {params.outdir}/mb2_abund_list.txt
                else
                    echo {params.outdir}/mb2_${sample%.*}.txt >> {params.outdir}/mb2_abund_list.txt
                fi
            done

            mkdir -p {params.outdir}/maxbin2_out/

            # Run maxbin2
            run_MaxBin.pl \
            -contig {input.contigs} \
            -markerset 107 \
            -thread {params.threads} \
            -min_contig_length 1500 \
	        -out {params.outdir}/maxbin2_out/bin \
        	-abund_list {params.outdir}/mb2_abund_list.txt


            # Create output for the next rule
            touch {output}

        fi
        """