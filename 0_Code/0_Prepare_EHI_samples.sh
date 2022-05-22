################################################################################
################################################################################
################################################################################
# This BASH script prepares the 1_Preprocess_QC.snakefile for EHI sample input.
# Raphael Eisenhofer 4/2022
#         .----------------.  .----------------.  .----------------.
#        | .--------------. || .--------------. || .--------------. |
#        | |  _________   | || |  ____  ____  | || |     _____    | |
#        | | |_   ___  |  | || | |_   ||   _| | || |    |_   _|   | |
#        | |   | |_  \_|  | || |   | |__| |   | || |      | |     | |
#        | |   |  _|  _   | || |   |  __  |   | || |      | |     | |
#        | |  _| |___/ |  | || |  _| |  | |_  | || |     _| |_    | |
#        | | |_________|  | || | |____||____| | || |    |_____|   | |
#        | |              | || |              | || |              | |
#        | '--------------' || '--------------' || '--------------' |
#         '----------------'  '----------------'  '----------------'
################################################################################
################################################################################
################################################################################

## Run preprocessing pipeline
#Rename raw sequencing files to EHIXXXXX numbers
#Create read groups (extract species column, sort, get unique values)
cut -f5 SEB001_experiment_checklist.tsv | sed '1d; s/ /_/g' | sort | uniq > genome_groups.tsv
while read group;
  do mkdir -p /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/$group &&
     mkdir -p /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/1_References/$group;
   done < genome_groups.tsv

#Manually put reference genomes in their respective folders
#Can automate this eventually, as AirTable has ftp links (or with preindexed computerome paths)

#extract first (EHIXXXXX#), fifth (species), and last (reverse file name)
cut -f1,5,15 SEB001_experiment_checklist.tsv | sed '1d; s/ /_/g' > EHIno_group_filename.tsv

#Sym link raw reads
ln -s /home/projects/ku-cbd/people/antalb/ehi_novogene/* /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed

#Create folders in 1_Untrimmed, rename sym-linked files
while read EHI group filename;
  do mv /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/*"$filename" /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/"$group"/"$EHI"_1.fastq.gz &&
     mv /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/*${filename/_1.fq/_2.fq} /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/"$group"/"$EHI"_2.fastq.gz;
done < EHIno_group_filename.tsv



##Setup a snakefile per species
mkdir temp_snakefiles

for group in /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/1_References/*;
  do cp /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/0_Code/1_Preprocess_QC.snakefile temp_snakefiles/$(basename $group)_1_Preprocess_QC.snakefile;
done

# Edit sample group paths for read input
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@1_Untrimmed/\*_1.fastq.gz@1_Untrimmed/$(basename ${group/_1_Preprocess_QC.snakefile/})/*_1.fastq.gz@" $group;
done

for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@1_Untrimmed/{sample}_1.fastq.gz@1_Untrimmed/$(basename ${group/_1_Preprocess_QC.snakefile/})/{sample}_1.fastq.gz@" $group;
done

for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@1_Untrimmed/{sample}_2.fastq.gz@1_Untrimmed/$(basename ${group/_1_Preprocess_QC.snakefile/})/{sample}_2.fastq.gz@" $group;
done


# Edit target reference genome path for input/output in rule 'index_ref' and input in rule 'map_to_ref'
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@\"1_References@\"1_References/$(basename ${group/_1_Preprocess_QC.snakefile/})@" $group;
done

# Edit snakefile to create a 'log' and 'benchmark' folder per group
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@0_Logs@0_Logs/$(basename ${group/_1_Preprocess_QC.snakefile/})@" $group;
done


##Launch snakefiles on Computerome!
for i in temp_snakefiles/*.snakefile;
  do snakemake -s snakefile -j 30 --cluster "qsub -l nodes=1:ppn={threads},{resources.mem_gb}G,walltime=00:NN:00:00" --use-conda --conda-frontend conda;
done
