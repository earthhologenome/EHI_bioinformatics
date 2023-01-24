#!/bin/bash

################################################################################
################################################################################
################################################################################
# This BASH script uploads preprocessed EHI reads and BAMs to the ENA.
# Raphael Eisenhofer 1/2023
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


# Helper function
help() {
    echo "###############################################"
    echo "###############################################"
    echo "Usage: $0 [-m metadata] [-i input_files] [-o output_xmls] [-c credentials]"
    echo " -m metadata = your metadata file, e.g. path/to/metadata.tsv"
    echo " -i input_files = path to your analysis files, e.g. path/to/files/"
    echo " -o ouput_xmls = path to where you wish to save .xmls and recipts, e.g. path/to/output/"
    echo " -u username = your ENA username, e.g. Webin-13337"
    echo " -p password = your ENA password, e.g. aWeSoMePaSsWoRd1!"
}

# Load in variables
while getopts ":m:i:o:c:u:p" opt; do
    case $opt in
        m) metadata="$OPTARG";;
        i) input_files="$OPTARG";;
        o) output_xmls="$OPTARG";;
        u) username="$OPTARG";;
        p) password="$OPTARG";;
        \?) help; exit 1;;
        :) echo "Option -$OPTARG requires an argument."; exit 1;;
    esac
done

# If no input -> show help message and exit
if [ -z "$metadata" ] || [ -z "$input_files" ] || [ -z "$output_xmls" ] || [ -z "$username" ] || [ -z "$password" ]; then
    help
    exit 1
fi


# Upload the analysis files (.bams & .fqs) to the ENA holding zone:
echo "Uploading analysis files to the ENA data holding zone, please wait..."

# first we have to set an environmental variable with our password so ascp does not prompt us:
export ASPERA_SCP_PASS=$password

# n.b. ' | tee aspera_log.txt' saves stdout to file for troubleshooting, while still printing it to stdout for user feedback.
ascp -QT -l300M -L- `pwd`"$input_files"/* "$username":@webin.ebi.ac.uk: | tee aspera_log.txt

echo "DONE!"

# loop through each line of the metadata
echo "Creating analysis XML files..."

# create a hidden backup metadata table, remove header of metadata, as it isn't required
cp $metadata ".$metadata_BACKUP"
sed -i'' '1d' $metadata

mkdir -p $output_xmls

while IFS=$'\t' read -r alias analysis_type study_ref sample_ref experiment_ref run_ref analysis_code reference_genome project analysis_center analysis_date file_name1 file_name2 file_type checksum_method checksum1 checksum2 assay_type analysis_protocol ; do

# if statement for file type (BAM / fastq)
if [ $file_type == "bam" ]
then


  # create the xml file for each sample (BAM)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name1}.xml"
  echo "    <ANALYSIS alias="\"$alias\"">" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <FILE filename="\"$file_name1\"" filetype="\"$file_type\""" >> "$output_xmls"/"${file_name1}.xml"
  echo "                checksum_method="\"$checksum_method\"" checksum="\"$checksum1\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$assay_type</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name1}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name1}.xml"


else


  # create the xml file for each sample (fastq)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name1}.xml"
  echo "    <ANALYSIS alias="\"$alias\"">" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <FILE filename="\"$file_name1\"" filetype="\"$file_type\""" >> "$output_xmls"/"${file_name1}.xml"
  echo "                checksum_method="\"$checksum_method\"" checksum="\"$checksum1\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <FILE filename="\"$file_name2\"" filetype="\"$file_type\""" >> "$output_xmls"/"${file_name1}.xml"
  echo "                checksum_method="\"$checksum_method\"" checksum="\"$checksum2\""/>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$assay_type</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name1}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name1}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name1}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name1}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name1}.xml"


fi

done < "$metadata"

echo "DONE!"


# Create MD5 hashes of files, then loop through each XML replacing with the hash for each file/s
echo "Creating MD5 hashes for analysis files and updating .XMLs..."

# create md5s:
for i in "$input_files"/*.bam;
    do md5sum $i > ${i/.bam/.md5};
done

for i in "$input_files"/*.fq.gz;
    do md5sum $i > ${i/.fq.gz/.md5};
done


# replace placeholder variables in XMLs with actual md5 hashes:
for i in "$input_files"/*.bam;
    do xml="$output_xmls/$(basename "$i".xml)"
       hash1=$(cut -f1 -d ' ' ${i/.bam/.md5})
       sed -i "s/PLACEHOLDER1/$hash1/g" $xml;
done

for i in "$input_files"/*_1.fq.gz;
    do xml="$output_xmls/$(basename "$i".xml)"
       hash1=$(cut -f1 -d ' ' ${i/.fq.gz/.md5})
       hash2=$(cut -f1 -d ' ' ${i/_1.fq.gz/_2.md5})
       sed -i "s/PLACEHOLDER1/$hash1/g" $xml
       sed -i "s/PLACEHOLDER2/$hash2/g" $xml;
done

echo "DONE!"

# Create analysis submission XML (only once)
echo "Creating submission XML file..."

echo "<SUBMISSION>" >> submission.xml
echo "   <ACTIONS>" >> submission.xml
echo "      <ACTION>" >> submission.xml
echo "         <ADD/>" >> submission.xml
echo "      </ACTION>" >> submission.xml
echo "   </ACTIONS>" >> submission.xml
echo "</SUBMISSION>" >> submission.xml


# Loop over each analysis XML, submitting them to the ENA and saving the receipt
echo "Submitting XML files to the ENA..."

for i in $output_xmls/*.xml;
    do curl -u "$username":"$password" -F "SUBMISSION=@submission.xml" -F "ANALYSIS=@$i" "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/";
done > $output_xmls/${i/.xml/_RECEIPT.xml}


echo "DONE! Have a good day :-)"
echo " "
echo "Your XML files and submission receipts are here:"
echo `pwd`/$output_xmls
