####
# Split count and coverage
# By Antton Alberdi (antton.alberdi@sund.ku.dk)
# 16/06/2023
####

#Load libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

#Parse input-output
parser <- ArgumentParser()
parser$add_argument("-i", "--input", action="store", help="Input combined file")
parser$add_argument("-c", "--counts", action="store", help="Count table file")
parser$add_argument("-v", "--coverage", action="store", help="Coverage table file")
args <- parser$parse_args()

# Split table

inputtable <- read.table(args$input, header = TRUE, sep = "\t")

inputtable %>%
    select(Genome, contains("Read.Count")) %>%
    rename_with(~ str_remove(., ".Read.Count"), contains("Read.Count")) %>%
    write.table(., file=args$counts, sep = "\t", quote = FALSE, row.names = FALSE)

inputtable %>%
    select(Genome, contains("Covered.Fraction")) %>%
    rename_with(~ str_remove(., ".Covered.Fraction"), contains("Covered.Fraction")) %>%
    write.table(., file=args$coverage, sep = "\t", quote = FALSE, row.names = FALSE)