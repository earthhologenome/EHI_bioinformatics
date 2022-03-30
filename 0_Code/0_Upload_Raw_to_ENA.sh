################################################################################
################################################################################
################################################################################
# This BASH script automatically uploads raw EHI reads to the ENA.
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

### General setup
# Create temp dir, get reference genome groups, create sample group directories
mkdir temp
cut -f3,4 -d ',' PR_Preprocessing_Input.csv | sed '1d;' | uniq | tr ' ' '_' | tr ',' '\t' > temp/group_genome.txt
cut -f1,3,10,11 -d ',' PR_Preprocessing_Input.csv | sed '1d;' | tr ' ' '_' | tr ',' '\t' | sed 's/_1.fq.gz//g' > temp/group_alias_file_seqbatch.txt

while read group ref;
  do mkdir -p 2_Reads/1_Untrimmed/$group &&
     mkdir -p 1_References/$group;
   done < temp/group_genome.txt

# Download reference genomes, and put them into their respective sample group folders
## NOTE -> ADD IF STATEMENT TO NOT RUN IF REF FILE ALREADY EXISTS?
while reads group ref;
  do wget $ref &&
     mv $ref 1_References/$group/;
   done < temp/group_genome.txt

# Rename raw filenames to EHI alias (e.g. EHI00001)
## NOTE -> Option to add hanlding of disparate file suffixes (e.g. .fq and .fastq)
while read alias group file seqbatch;
  do cp 2_Reads/0_Raw/$seqbatch/"$file"_1.fq.gz 2_Reads/1_Untrimmed/$group/"$alias"_R_1.fastq.gz &&
     cp 2_Reads/0_Raw/$seqbatch/"$file"_2.fq.gz 2_Reads/1_Untrimmed/$group/"$alias"_R_2.fastq.gz;
   done < temp/group_alias_file_seqbatch.txt


### Setup 1_Preprocess_QC snakefiles -- one per sample group
for group in 1_References/*;
  do cp 0_Code/1_Preprocess_QC.snakefile temp/$(basename $group)_1_Preprocess_QC.snakefile;
   done

# Edit sample group paths for read input
for group in temp/*.snakefile;
  do sed -i'' "s@1_Untrimmed/\*_1.fq.gz@1_Untrimmed/$(basename ${group/_1_Preprocess_QC.snakefile/})/*_1.fq.gz@" $group;
   done

# Edit target reference genome path for input/output in rule 'index_ref' and input in rule 'map_to_ref'
for group in temp/*.snakefile;
  do sed -i'' "s@\"1_References@\"1_References/$(basename ${group/_1_Preprocess_QC.snakefile/})@" $group;
   done

# Edit snakefile to create a 'log' and 'benchmark' folder per group
for group in temp/*.snakefile;
  do sed -i'' "s@0_Logs@0_Logs/$(basename ${group/_1_Preprocess_QC.snakefile/})@" $group;
   done


# Clean up temp files
rm -r temp
