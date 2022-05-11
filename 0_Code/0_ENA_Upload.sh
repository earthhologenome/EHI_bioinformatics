## First time, setting up conda environment and installing forked ena-upload-cli version:
git clone https://github.com/EisenRa/ena-upload-cli.git

#Create conda environment:
conda env create --prefix .ENAup --file ena_upload_conda.yml

#Install from source
cd ena-upload-cli
python setup.py install

## Setup ENA studies first, if needed


## Setup ENA specimen accessions
#Rename file by sequencing batch
mv Captures-ENA_template_specimens_ERC000053.csv SEB001_Captures-ENA_template_specimens_ERC000053.csv

#Convert csv to tsv
cat SEB001_Captures-ENA_template_specimens_ERC000053.csv | tr ',' '\t' > SEB001_Captures-ENA_template_specimens_ERC000053.tsv

#Submit to ENA
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample SEB001_Captures-ENA_template_specimens_ERC000053.tsv \
--secret /home/projects/ku-cbd/people/rapeis/EHI/0_Software/.secret.yml \
--checklist ERC000053

## Grab ENA generated SAMPLE ACCESSIONS (SAME######)
#This needs to be joined into the sample sheets for ENA compatibility
mv receipt.xml SEB001_receipt_specimen.xml
grep 'alias=' SEB001_receipt_specimen.xml | cut -d '"' -f4 > SEB001_ENA_specimen_ALIASES.tsv
grep 'SAMEA' SEB001_receipt_specimen.xml | cut -d '"' -f2 > SEB001_ENA_specimen_SAMEA_ACCESSIONS.tsv
paste SEB001_ENA_ALIASES.tsv SEB001_ENA_SAMEA_ACCESSIONS.tsv > SEB001_ENA_EHI_specimen_mapping.tsv

#Next, download the mapping file and add the ENA sample accession to the air table
#Make sure you sort both the AirTable and tsv to ensure proper mapping. DOUBLE CHECK!
scp co:/home/projects/ku-cbd/people/rapeis/EHI/SEB001/ENA_upload/SEB001_ENA_EHI_mapping.tsv .

#Once you've done that, export the next sample accession sheets from AirTable and
#upload to the working directory.
mv Samples-ENA_template_samples_ERC000013.csv SEB001_Samples-ENA_template_samples_ERC000013.csv
mv Samples-ENA_template_samples_ERC000053.csv SEB001_Samples-ENA_template_samples_ERC000053.csv
cat SEB001_Samples-ENA_template_samples_ERC000013.csv | tr ',' '\t' > SEB001_Samples-ENA_template_samples_ERC000013.tsv
cat SEB001_Samples-ENA_template_samples_ERC000053.csv | tr ',' '\t' > SEB001_Samples-ENA_template_samples_ERC000053.tsv

## Setup ENA sample accessions (gut metagenome ERC000013 and/or host genome ERC000053)
#Host genome ERC000053
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample SEB001_Samples-ENA_template_samples_ERC000013.tsv \
--secret /home/projects/ku-cbd/people/rapeis/EHI/0_Software/.secret.yml \
--checklist ERC000013

mv receipt.xml SEB001_receipt_sample_ERC000013.xml
grep 'alias=' SEB001_receipt_sample_ERC000013.xml | cut -d '"' -f4 > SEB001_ENA_ERC000013_ALIASES.tsv
grep 'SAMEA' SEB001_receipt_sample_ERC000013.xml | cut -d '"' -f2 > SEB001_ENA_ERC000013_SAMEA_ACCESSIONS.tsv
paste SEB001_ENA_ERC000013_ALIASES.tsv SEB001_ENA_ERC000013_SAMEA_ACCESSIONS.tsv > SEB001_ENA_EHI_ERC000013_mapping.tsv

#Host associated metagenome ERC000013
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--sample SEB001_Samples-ENA_template_samples_ERC000053.tsv \
--secret /home/projects/ku-cbd/people/rapeis/EHI/0_Software/.secret.yml \
--checklist ERC000053

mv receipt.xml SEB001_receipt_sample_ERC000053.xml
grep 'alias=' SEB001_receipt_sample_ERC000053.xml | cut -d '"' -f4 > SEB001_ENA_ERC000053_ALIASES.tsv
grep 'SAMEA' SEB001_receipt_sample_ERC000053.xml | cut -d '"' -f2 > SEB001_ENA_ERC000053_SAMEA_ACCESSIONS.tsv
paste SEB001_ENA_ERC000053_ALIASES.tsv SEB001_ENA_ERC000053_SAMEA_ACCESSIONS.tsv > SEB001_ENA_EHI_ERC000053_mapping.tsv

#Now make sure to add these sample accesions to the AirTable:
#Make sure you sort both the AirTable and tsv to ensure proper mapping. DOUBLE CHECK!
scp co:/home/projects/ku-cbd/people/rapeis/EHI/SEB001/ENA_upload/SEB001_ENA_EHI_ERC* .


## Setup ENA experiment and run accessions, upload raw reads:
#Export the 'ENA_experiment_checklist' from the EHI AirTable, upload to working directory
mv SE\ \(Samples\)-ENA_experiment_checklist.csv SEB001_experiment_checklist.csv
cat SEB001_experiment_checklist.csv | tr ',' '\t' > SEB001_experiment_checklist.tsv

###########################################################
#Run the first line of code in the 0_Prepare_EHI_samples.sh
###########################################################

##Create the ENA run template file:
#Create header
echo -e "alias\texperiment_alias\tfile_name\tfile_type\tfile_checksum" > run_headers.tsv
#Pull out EHI number
cut -f1 SEB001_experiment_checklist.tsv | sed '1d;' > EHInumbers.tsv
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
cat run_headers.tsv merged.tsv > SEB001_ENA_run_sheet.tsv


##Submit the experiment and run ENA sheets, and upload the raw reads to ENA!
ena-upload-cli \
--action add \
--center 'Earth Hologenome Initiative' \
--experiment SEB001_experiment_checklist.tsv \
--run SEB001_ENA_run_sheet.tsv \
--secret /home/projects/ku-cbd/people/rapeis/EHI/0_Software/.secret.yml \
--data /home/projects/ku-cbd/people/rapeis/EHI/SEB001/EHI_bioinformatics/2_Reads/1_Untrimmed/*/*.fastq.gz \
--draft
