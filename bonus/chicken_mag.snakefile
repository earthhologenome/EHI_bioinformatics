################################################################################
################################################################################
################################################################################
### Snakefile for updating chicken MAG catalogues with new MAGs
### Raphael Eisenhofer 04/2023

## Load dependencies, set MAG variable
import os
import glob

MAG = [os.path.basename(fn).replace(".fa.gz", "")
            for fn in glob.glob(os.path.join(f"new_mags/", "*.fa.gz"))]

print("Detected the following new MAGs:")
print(MAG)

## Declare our desired output
rule all:
    input:
        "upload_complete"

###############################################################################
## Download dereplicated MAGs from ERDA
rule get_derepd_mags:
    output:
        mags = directory("dereplicated_mags/")
        mag_stats = "mag_stats.csv"
    conda:
    threads: 1
    shell:
        """
        lftp sftp://erda -e "mirror -R dereplicated_mags/ /chicken_mag_catalogue/dereplicated/; bye"
        lftp sftp://erda -e "get /chicken_mag_catalogue/mag_stats.csv; bye"
        """

rule checkm_mags:
    output:
        new_mag_stats = "new_mag_stats.csv"
    conda:
    threads: 16
    shell:
        """
        checkm 
        """    

rule dereplicate:
    input:
        "new_mag_stats.csv"
    output:
        "dRep/mag_stats.csv"
    conda:
    threads: 16
    shell:
        """
        # move new mags into previously dereplicated folder
        mv new_mags/*.fa.gz dereplicated_mags/

        # run dRep
        dRep dereplicate dRep \
            -p {threads} \
            -sa 0.98 \
            -comp 50 \
            -con 10 \
            -g mags/*.fa.gz

        # create new 'mag_stats.csv' for upload to erda
        
        """    

rule upload_to_erda:
    input:
        "dRep/mag_stats.csv"
    output:
        "upload_complete"
    conda:

    threads: 1
    shell:
        """
        # remove previous dereplicated genomes and mag_stats.csv from erda

        
        """