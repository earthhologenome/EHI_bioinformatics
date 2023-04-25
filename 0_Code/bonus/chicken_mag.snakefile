################################################################################
################################################################################
################################################################################
### Snakefile for updating chicken MAG catalogues with new MAGs
### Raphael Eisenhofer 04/2023

## Load dependencies, set MAG variable
import os
import glob

MAG = [os.path.basename(fn).replace(".fa.gz", "")
            for fn in glob.glob(os.path.join(f"mags/", "*.fa.gz"))]

print("Detected the following MAGs:")
print(MAG)

## Declare our desired output
rule all:
    input:

###############################################################################
## Download dereplicated MAGs from ERDA
rule get_derepd_mags:
    input:
    output:
        directory("dereplicated_mags/")
    shell:
        """

    
        """