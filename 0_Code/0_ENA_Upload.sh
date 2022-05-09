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
echo "sample derived from" > SEB001_ENA_SAMEA_ACCESSIONS.tsv
grep 'SAMEA' SEB001_receipt_specimen.xml | cut -d '"' -f2 >> SEB001_ENA_SAMEA_ACCESSIONS.tsv

join
paste ENA_template_samples_ERC000053.tsv ENA_SAMEA_ACCESSIONS.tsv > ENA_template_samples_ERC000053_SAMEA.tsv
paste ENA_template_samples_ERC000013.tsv ENA_SAMEA_ACCESSIONS.tsv > ENA_template_samples_ERC000013_SAMEA.tsv

## Setup ENA sample accessions (gut metagenome ERC000013 and/or host genome ERC000053)


mv receipt.xml receipt_experiment.xml
