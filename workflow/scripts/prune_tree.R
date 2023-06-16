####
# Prune phylogenetic tree for EHI DM output
# By Antton Alberdi (antton.alberdi@sund.ku.dk)
# 16/06/2023
####

#Load libraries
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(argparse))

#Parse input-output
parser <- ArgumentParser()
parser$add_argument("-i", "--input", action="store", help="Input tree file")
parser$add_argument("-o", "--output", action="store", help="Output tree file")
args <- parser$parse_args()

#Prune and output tree
input_tree <- ape::read.tree(args$input)
input_tips <- input_tree$tip.label
output_tips <- input_tips[grepl("^EHA", input_tips)]
output_tree <- ape::keep.tip(input_tree,output_tips)
ape::write.tree(output_tree,args$output)