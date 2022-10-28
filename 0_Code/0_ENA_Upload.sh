################################################################################
################################################################################
################################################################################
# This BASH script uploads raw EHI sequencing reads and metadata to the ENA.
# Raphael Eisenhofer 5/2022
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

# Setup variables for the pipeline

echo "Please type in the sequencing batch number -- e.g. SEB001"
read SeqBatch

echo "Are you adding new host genomic samples? 'yes' or 'no'"
read new_genomic_sample

## First time, setting up conda environment and installing forked ena-upload-cli version:
#git clone https://github.com/EisenRa/ena-upload-cli.git

#Create conda environment:
#conda env create --prefix .ENAup --file ena_upload_conda.yml

#Install from source
#cd ena-upload-cli
#python setup.py install

## Setup ENA studies first, if needed

## Setup ENA specimen accessions
#Rename file by sequencing batch
mv Captures-ENA_template_specimens_ERC000053.csv "$SeqBatch"_Captures-ENA_template_specimens_ERC000053.csv

#Convert csv to tsv
cat "$SeqBatch"_Captures-ENA_template_specimens_ERC000053.csv | tr ',' '\t' > "$SeqBatch"_Captures-ENA_template_specimens_ERC000053.tsv

#Submit to ENA
#N.B. the .secret.yml is currently stored in Raphael's personal folder
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample "$SeqBatch"_Captures-ENA_template_specimens_ERC000053.tsv \
--secret /projects/mjolnir1/people/ncl550/0_software/.secret.yml \
--checklist ERC000053

## Grab ENA generated SAMPLE ACCESSIONS (SAME######)
#This needs to be joined into the sample sheets for ENA compatibility
mv receipt.xml "$SeqBatch"_receipt_specimen.xml
grep 'alias=' "$SeqBatch"_receipt_specimen.xml | cut -d '"' -f4 > "$SeqBatch"_ENA_specimen_ALIASES.tsv
grep 'SAMEA' "$SeqBatch"_receipt_specimen.xml | cut -d '"' -f2 > "$SeqBatch"_ENA_specimen_SAMEA_ACCESSIONS.tsv
paste "$SeqBatch"_ENA_ALIASES.tsv "$SeqBatch"_ENA_SAMEA_ACCESSIONS.tsv > "$SeqBatch"_ENA_EHI_specimen_mapping.tsv


## This is for the host genomic samples
if [ $new_genomic_sample == "yes" ]
then

mv Samples-ENA_template_samples_ERC000053.csv "$SeqBatch"_Samples-ENA_template_samples_ERC000053.csv
cat "$SeqBatch"_Samples-ENA_template_samples_ERC000053.csv | tr ',' '\t' > "$SeqBatch"_Samples-ENA_template_samples_ERC000053.tsv

## Setup ENA sample accessions (gut metagenome ERC000013 and/or host genome ERC000053)

#Host genome ERC000053
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample "$SeqBatch"_Samples-ENA_template_samples_ERC000053.tsv \
--secret /projects/mjolnir1/people/ncl550/0_software/.secret.yml \
--checklist ERC000053

mv receipt.xml "$SeqBatch"_receipt_sample_ERC000053.xml
grep 'alias=' "$SeqBatch"_receipt_sample_ERC000053.xml | cut -d '"' -f4 > "$SeqBatch"_ENA_ERC000053_ALIASES.tsv
grep 'SAMEA' "$SeqBatch"_receipt_sample_ERC000053.xml | cut -d '"' -f2 > "$SeqBatch"_ENA_ERC000053_SAMEA_ACCESSIONS.tsv
paste "$SeqBatch"_ENA_ERC000053_ALIASES.tsv "$SeqBatch"_ENA_ERC000053_SAMEA_ACCESSIONS.tsv > "$SeqBatch"_ENA_EHI_ERC000053_mapping.tsv

else
echo "No genomic samples are being added."

fi

#Host associated metagenome ERC000013
mv Samples-ENA_template_samples_ERC000013.csv "$SeqBatch"_Samples-ENA_template_samples_ERC000013.csv
cat "$SeqBatch"_Samples-ENA_template_samples_ERC000013.csv | tr ',' '\t' > "$SeqBatch"_Samples-ENA_template_samples_ERC000013.tsv

ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample "$SeqBatch"_Samples-ENA_template_samples_ERC000013.tsv \
--secret /projects/mjolnir1/people/ncl550/0_software/.secret.yml \
--checklist ERC000013

mv receipt.xml "$SeqBatch"_receipt_sample_ERC000013.xml
grep 'alias=' "$SeqBatch"_receipt_sample_ERC000013.xml | cut -d '"' -f4 > "$SeqBatch"_ENA_ERC000013_ALIASES.tsv
grep 'SAMEA' "$SeqBatch"_receipt_sample_ERC000013.xml | cut -d '"' -f2 > "$SeqBatch"_ENA_ERC000013_SAMEA_ACCESSIONS.tsv
paste "$SeqBatch"_ENA_ERC000013_ALIASES.tsv "$SeqBatch"_ENA_ERC000013_SAMEA_ACCESSIONS.tsv > "$SeqBatch"_ENA_EHI_ERC000013_mapping.tsv


## Setup ENA experiment and run accessions, upload raw reads:
#Export the 'ENA_experiment_checklist' from the EHI AirTable, upload to working directory
mv SE\ \(Samples\)-ENA_experiment_checklist.csv "$SeqBatch"_experiment_checklist.csv
cat "$SeqBatch"_experiment_checklist.csv | tr ',' '\t' > "$SeqBatch"_experiment_checklist.tsv


#Rename raw sequencing files to EHIXXXXX numbers
#Create read groups (extract species column, sort, get unique values)
cut -f5 "SeqBatch"_experiment_checklist.tsv | sed '1d; s/ /_/g' | sort | uniq > genome_groups.tsv
while read group;
  do mkdir -p `pwd`/2_Reads/1_Untrimmed/$group &&
     mkdir -p `pwd`/1_References/$group;
   done < genome_groups.tsv

#Manually put reference genomes in their respective folders
#Can automate this eventually, as AirTable has ftp links (or with preindexed computerome paths)

#extract first (EHIXXXXX#), fifth (species), and last (reverse file name)
cut -f1,5,15 "SeqBatch"_experiment_checklist.tsv | sed '1d; s/ /_/g' > EHIno_group_filename.tsv

#Move the read files into their appropriate folders
ln -s `pwd`/2_Reads/0_Raw/* `pwd`/2_Reads/1_Untrimmed

#Create folders in 1_Untrimmed, rename sym-linked files
while read EHI group filename;
  do mv `pwd`/2_Reads/1_Untrimmed/*"$filename" `pwd`/2_Reads/1_Untrimmed/"$group"/"$EHI"_1.fastq.gz &&
     mv `pwd`/2_Reads/1_Untrimmed/*${filename/_1.fq/_2.fq} `pwd`/2_Reads/1_Untrimmed/"$group"/"$EHI"_2.fastq.gz;
done < EHIno_group_filename.tsv


##Setup a snakefile per species
mkdir -p temp_snakefiles

for group in `pwd`/1_References/*;
  do cp `pwd`/0_Code/1_Preprocess_QC.snakefile temp_snakefiles/$(basename $group)_1_Preprocess_QC.snakefile;
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

# Edit snakefile to have group name in final report
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@_report.tsv@$(basename ${group/_1_Preprocess_QC.snakefile/})_report.tsv@g" $group;
done

# Edit snakefile to have group name in nonpareil metadata
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@_nonpareil_metadata@$(basename ${group/_1_Preprocess_QC.snakefile/})_nonpareil_metadata@g" $group;
done

# Edit snakefile to output non-host reads into group-specific folders
for group in temp_snakefiles/*.snakefile;
  do sed -i'' "s@2_Reads/4_Host_removed/@2_Reads/4_Host_removed/$(basename ${group/_1_Preprocess_QC.snakefile/})/@g" $group;
done




##Create the ENA run template file:
#Create header
echo -e "alias\texperiment_alias\tfile_name\tfile_type\tfile_checksum" > run_headers.tsv
#Pull out EHI number
cut -f1 "$SeqBatch"_experiment_checklist.tsv | sed '1d;' > EHInumbers.tsv
#Create a copy (we need 2 per EHI number -- R1/R2)
cat EHInumbers.tsv EHInumbers.tsv > EHInumbers2.tsv
#Add _1 suffix to EHI number???
#Add raw R1 fastq filename suffix to EHI number
sed 's/$/_1.fastq.gz/g' EHInumbers.tsv > R1.tsv
#Add raw R2 fastq filename suffix to EHI number
sed 's/$/_2.fastq.gz/g' EHInumbers.tsv > R2.tsv
#Create column containing 'fastq' string
while read EHI; do echo "fastq" >> fastq.tsv; done < EHInumbers2.tsv
#Merge table
cat R1.tsv R2.tsv > R1R2.tsv
paste EHInumbers2.tsv EHInumbers2.tsv R1R2.tsv fastq.tsv > merged.tsv
cat run_headers.tsv merged.tsv > "$SeqBatch"_ENA_run_sheet.tsv

# ##In this example, there is no EHI00048, so remove from exp/run table
# cp "$SeqBatch"_experiment_checklist.tsv "$SeqBatch"_experiment_checklist.tsv1 && rm "$SeqBatch"_experiment_checklist.tsv
# grep -v 'EHI00048' "$SeqBatch"_experiment_checklist.tsv1 > "$SeqBatch"_experiment_checklist.tsv && rm "$SeqBatch"_experiment_checklist.tsv1
# cp "$SeqBatch"_ENA_run_sheet.tsv "$SeqBatch"_ENA_run_sheet.tsv1 && rm "$SeqBatch"_ENA_run_sheet.tsv
# grep -v 'EHI00048' "$SeqBatch"_ENA_run_sheet.tsv1 > "$SeqBatch"_ENA_run_sheet.tsv && rm "$SeqBatch"_ENA_run_sheet.tsv1


##FTP is not very stable (had a lot of issues with stalling halfway through submission)
##We'll use aspera cli instead:
# https://d3gcli72yxqn2z.cloudfront.net/connect_latest/v4/bin/ibm-aspera-connect_4.1.3.93_linux.tar.gz
# or: https://www.ibm.com/aspera/connect/

##Upload the files via aspera (they are stored temporarily on ENA server, and can be linked after)
ascp -QT -l100M -L- `pwd`/2_Reads/1_Untrimmed/*/*.fastq.gz Webin-61693@webin.ebi.ac.uk:

##Submit the experiment and run ENA sheets, linking the files we uploaded in the previous step
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--experiment "$SeqBatch"_experiment_checklist.tsv \
--run "$SeqBatch"_ENA_run_sheet.tsv \
--secret /projects/mjolnir1/people/ncl550/0_software/.secret.yml \
--data `pwd`/2_Reads/1_Untrimmed/*/*.fastq.gz \
--no_data_upload

mv receipt.xml "$SeqBatch"_receipt_run.xml

##Now manually enter in the ERX & ERR accessions into the EHI AirTable
grep '<RUN' "$SeqBatch"_receipt_run.xml | cut -f4 -d '"' > "$SeqBatch"_run_aliases.tsv
grep '<EXPERIMENT' "$SeqBatch"_receipt_run.xml | cut -f2 -d '"' > "$SeqBatch"_run_ERX.tsv
grep '<RUN' "$SeqBatch"_receipt_run.xml | cut -f2 -d '"' > "$SeqBatch"_run_ERR.tsv
paste "$SeqBatch"_run_aliases.tsv "$SeqBatch"_run_ERR.tsv "$SeqBatch"_run_ERX.tsv > "$SeqBatch"_run_mapping.tsv
