####
# Prune phylogenetic tree for EHI DM output
# By Antton Alberdi (antton.alberdi@sund.ku.dk)
# 30/06/2023
# bug fix by Raphael 5/7/2023
####

#Load libraries
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


#Parse input-output
parser <- ArgumentParser()
parser$add_argument("-a", "--archaea", action="store", help="Input archaea tree file")
parser$add_argument("-b", "--bacteria", action="store", help="Input bacteria tree file")
parser$add_argument("-i", "--info", action="store", help="Genome info table")
parser$add_argument("-o", "--output", action="store", help="Output tree file")
args <- parser$parse_args()

#Load and bind trees
archaea_tree <- ape::read.tree(args$archaea)
bacteria_tree <- ape::read.tree(args$bacteria)
prokaryote_tree <- bind.tree(archaea_tree, bacteria_tree)

#Prune tree
genomes <- read.table(args$info, sep = '\t', skip = 1)[,1]
genomes <- str_replace_all(genomes, "$", ".fa")
output_tree <- ape::keep.tip(prokaryote_tree, genomes)

#Quality-check
if(length(genomes) == length(output_tree$tip.label)){
	ape::write.tree(output_tree, args$output)
}else{
	stop("Mag info and tree tips do not match")
}
