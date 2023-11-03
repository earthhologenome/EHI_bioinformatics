####
# Render final report
# By Antton Alberdi (antton.alberdi@sund.ku.dk)
# 16/06/2023
# Dependencies: pandoc
####

#Load libraries
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(R.utils))


#Parse input-output
parser <- ArgumentParser()
parser$add_argument("-s", "--script", action="store", help="Script file")
parser$add_argument("-b", "--batch", action="store", help="DM Batch code")
parser$add_argument("-p", "--pdf", action="store", help="PDF output file")
args <- parser$parse_args()

#Render documents
rmarkdown::render(args$script, params = list(batch=args$batch, count_file = paste0(args$batch,"_counts.tsv.gz"), coverage_file = paste0(args$batch,"_coverage.tsv.gz"), sample_file = paste0(args$batch,"_metadata.tsv.gz"), mags_file = paste0(args$batch,"_mag_info.tsv.gz"), tree_file = paste0(args$batch,".tree.gz"), kegg_file = paste0(args$batch,"_merged_kegg.tsv.gz")), output_format = "pdf_document", output_file = args$pdf)
#rmarkdown::render(args$script, params = list(count_file = "DMB0019_counts.tsv.gz", coverage_file = "DMB0019_coverage.tsv.gz", sample_file = "DMB0019_metadata.tsv.gz", mags_file = "DMB0019_mag_info.tsv.gz", tree_file = "DMB0019.tree.gz", kegg_file = paste0(args$batch,"_merged_kegg.tsv.gz"), output_format = "html_document", output_file = args$html)

#cd /Users/anttonalberdi/Downloads/DMB0043/
#Rscript final_report_render.r -s final_report.r -b DMB0043 -p DMB0043.pdf
