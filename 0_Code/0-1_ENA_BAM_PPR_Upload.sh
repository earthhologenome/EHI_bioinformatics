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
    echo " "
    echo "Usage: $0 [-m metadata] [-i input_files] [-o output_xmls] [-u username] [-p pass] [-t type]"
    echo " "
    echo " -m metadata = your metadata file, e.g. path/to/metadata.csv from AirTable"
    echo " -i input_files = path to your analysis files, e.g. path/to/files/"
    echo " -o ouput_xmls = path to where you wish to save .xmls and recipts, e.g. path/to/output/"
    echo " -u username = your ENA username, e.g. Webin-13337"
    echo " -p pass = your ENA password, e.g. aWeSoMePaSsWoRd1!"
    echo " -t type = what type of analysis, options= 'reads', 'bams', 'both'"
}

# Load in variables
while getopts ":m:i:o:u:p:t:" opt; do
    case $opt in
        m) metadata="$OPTARG";;
        i) input_files="$OPTARG";;
        o) output_xmls="$OPTARG";;
        u) username="$OPTARG";;
        p) pass="$OPTARG";;
        t) type="$OPTARG";;
        \?) help; exit 1;;
    esac
done

# If no input -> show help message and exit
if [ -z "$metadata" ] || [ -z "$input_files" ] || [ -z "$output_xmls" ] || [ -z "$username" ] || [ -z "$pass" ] || [ -z "$type" ]; then
    echo "All input fields are required!"
    help
    exit 1
fi

# Upload the analysis files (.bams & .fqs) to the ENA holding zone:
echo "Uploading analysis files to the ENA data holding zone, please wait..."

# first we have to set an environmental variable with our password so ascp does not prompt us:
export ASPERA_SCP_PASS=$pass

# n.b. ' | tee aspera_log.txt' saves stdout to file for troubleshooting, while still printing it to stdout for user feedback.
ascp -QT -l300M -L- `pwd`"$input_files"/* "$username":@webin.ebi.ac.uk: |& tee aspera_log.txt

echo "DONE!"

# loop through each line of the metadata
echo "Creating analysis XML files..."

# remove header of metadata, as it isn't required
sed '1d' $metadata > temp.tsv
    # stupid EOF issues!
    echo " " > eof.txt && cat temp.tsv eof.txt > $metadata
    rm temp.tsv && rm eof.txt

mkdir -p $output_xmls

# while IFS=$'\t' read -r PR_code alias_bam alias_fastq study_ref sample_ref experiment_ref run_ref analysis_code reference_genome project analysis_center analysis_date file_name_bam file_name_fastq1 file_name_fastq2 analysis_protocol ; do
while IFS=$',' read -r PR_code alias_bam alias_fastq study_ref sample_ref experiment_ref run_ref analysis_code reference_genome project analysis_center analysis_date file_name_bam file_name_fastq1 file_name_fastq2 analysis_protocol ; do

file_name_fastq=(${file_name_fastq1/_1.fq.gz})

# if type = bams, then only create bam XMLs
if [ $type == bams ]
then

  # create the xml file for each sample (BAM)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "    <ANALYSIS alias="\"$alias_bam\"">" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <FILE filename="\"$file_name_bam\"" filetype="\"bam\""" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM1\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>high throughput sequencing</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name_bam}.xml"


    # create md5s:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/*.bam;
            do md5sum $i | cut -f1 -d ' ' > ${i/.bam/.md5};
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/*.bam;
            do md5 $i | cut -f4 -d ' ' > ${i/.bam/.md5};
        done

    else
        echo "Operating system not supported!"
        exit 1
    fi

    # replace placeholder variables in XMLs with actual md5 hashes:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/*.bam;
            do xml="$output_xmls/$(basename "$i".xml)"
            hash1=$(cat ${i/.bam/.md5})
            sed -i "s/CHECKSUM1/$hash1/g" $xml;
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/*.bam;
            do xml="$output_xmls/$(basename "$i".xml)"
            hash1=$(cat ${i/.bam/.md5})
            sed -i '' "s/CHECKSUM1/$hash1/g" $xml;
        done
    fi


# IF type = reads, only create fastq XMLs
elif [ $type == reads ]
then

  # create the xml file for each sample (fastq)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "    <ANALYSIS alias="\"$alias_fastq\"">" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <FILE filename="\"$file_name_fastq1\"" filetype="\"fastq\""" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM2\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <FILE filename="\"$file_name_fastq2\"" filetype="\"fastq\""" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM3\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>high throughput sequencing</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name_fastq}.xml"

    # create md5s:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/*.fq.gz;
            do md5sum $i | cut -f1 -d ' ' > ${i/.fq.gz/.md5};
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/*.fq.gz;
            do md5 $i | cut -f4 -d ' ' > ${i/.fq.gz/.md5};
        done

    else
        echo "Operating system not supported!"
        exit 1
    fi

    # replace placeholder variables in XMLs with actual md5 hashes:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/"$file_name_fastq1";
            do xml="$output_xmls/$(basename "${i/_1.fq.gz/}".xml)"
            echo $file_name_fastq1
            echo $i
            echo $xml
            hash1=$(cat ${i/.fq.gz/.md5})
            hash2=$(cat ${i/_1.fq.gz/_2.md5})
            sed -i "s/CHECKSUM2/$hash1/g" $xml
            sed -i "s/CHECKSUM3/$hash2/g" $xml;
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/"$file_name_fastq1";
            do xml="$output_xmls/$(basename "${i/_1.fq.gz/}".xml)"
            hash1=$(cat ${i/.fq.gz/.md5})
            hash2=$(cat ${i/_1.fq.gz/_2.md5})
            sed -i '' "s/CHECKSUM2/$hash1/g" $xml
            sed -i '' "s/CHECKSUM3/$hash2/g" $xml;
        done
    fi


# IF type = both, create both types of XML file
elif [ $type == both ]
then

  # create the xml file for each sample (BAM)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "    <ANALYSIS alias="\"$alias_bam\"">" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <FILE filename="\"$file_name_bam\"" filetype="\"bam\""" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM1\""/>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>high throughput sequencing</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name_bam}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name_bam}.xml"


  # create the xml file for each sample (fastq)
  echo "<ANALYSIS_SET>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "    <ANALYSIS alias="\"$alias_fastq\"">" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <PROCESSED_READS/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </ANALYSIS_TYPE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <FILES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <FILE filename="\"$file_name_fastq1\"" filetype="\"fastq\""" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM2\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <FILE filename="\"$file_name_fastq2\"" filetype="\"fastq\""" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                checksum_method="\"MD5\"" checksum="\"CHECKSUM3\""/>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </FILES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Project</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$project</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>high throughput sequencing</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "    </ANALYSIS>" >> "$output_xmls"/"${file_name_fastq}.xml"
  echo "</ANALYSIS_SET>" >> "$output_xmls"/"${file_name_fastq}.xml"

    # create md5s:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/*.bam;
            do md5sum $i | cut -f1 -d ' ' > ${i/.bam/.md5};
        done

        for i in "$input_files"/*.fq.gz;
            do md5sum $i | cut -f1 -d ' ' > ${i/.fq.gz/.md5};
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/*.bam;
            do md5 $i | cut -f4 -d ' ' > ${i/.bam/.md5};
        done

        for i in "$input_files"/*.fq.gz;
            do md5 $i | cut -f4 -d ' ' > ${i/.fq.gz/.md5};
        done

    else
        echo "Operating system not supported!"
        exit 1
    fi


# replace placeholder variables in XMLs with actual md5 hashes:
    if [[ "$OSTYPE" == "linux"* ]]; then
        for i in "$input_files"/*.bam;
            do xml="$output_xmls/$(basename "$i".xml)"
            hash1=$(cat ${i/.bam/.md5})
            sed -i "s/CHECKSUM1/$hash1/g" $xml;
        done

        for i in "$input_files"/"$file_name_fastq1";
            do xml="$output_xmls/$(basename "${i/_1.fq.gz/}".xml)"
            hash1=$(cat ${i/.fq.gz/.md5})
            hash2=$(cat ${i/_1.fq.gz/_2.md5})
            sed -i "s/CHECKSUM2/$hash1/g" $xml
            sed -i "s/CHECKSUM3/$hash2/g" $xml;
        done

    elif [[ "$OSTYPE" == "darwin"* ]]; then
        for i in "$input_files"/*.bam;
            do xml="$output_xmls/$(basename "$i".xml)"
            hash1=$(cat ${i/.bam/.md5})
            sed -i '' "s/CHECKSUM1/$hash1/g" $xml;
        done

        for i in "$input_files"/"$file_name_fastq1";
            do xml="$output_xmls/$(basename "${i/_1.fq.gz/}".xml)"
            hash1=$(cat ${i/.fq.gz/.md5})
            hash2=$(cat ${i/_1.fq.gz/_2.md5})
            sed -i '' "s/CHECKSUM2/$hash1/g" $xml
            sed -i '' "s/CHECKSUM3/$hash2/g" $xml;
        done
    fi


fi

done < $metadata

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
    do curl -u "$username":"$pass" -F "SUBMISSION=@submission.xml" -F "ANALYSIS=@$i" "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/";
done > $output_xmls/${i/.xml/_RECEIPT.xml}


for i in 


success="false"

echo "DONE! Have a good day :-)"
echo " "
echo "Your XML files and submission receipts are here:"
echo `pwd`/$output_xmls
