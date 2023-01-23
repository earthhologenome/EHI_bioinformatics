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

# the input folder containing the reads or bams
input_files=$1
# the input metadata file for the samples (see X for template)
metadata=$2

# loop through each line of the metadata
while IFS=$'\t' read -r alias analysis_type study_ref sample_ref experiment_ref run_ref analysis_code reference_genome project analysis_center analysis_date file_names file_types checksum_method checksum assay_type analysis_protocol ; do
  # create the xml file for each sample
  echo "<ANALYSIS_SET>" >> "${alias}.xml"
  echo "    <ANALYSIS alias="\"$alias\"">" >> "${alias}.xml"
  echo "        <TITLE>Processed reads</TITLE>" >> "${alias}.xml"
  echo "        <DESCRIPTION>EHI preprocessed files</DESCRIPTION>" >> "${alias}.xml"
  echo "        <STUDY_REF accession="\"$study_ref\""/>" >> "${alias}.xml"
  echo "        <SAMPLE_REF accession="\"$sample_ref\""/>" >> "${alias}.xml"
  echo "        <EXPERIMENT_REF accession="\"$experiment_ref\""/>" >> "${alias}.xml"
  echo "        <RUN_REF accession="\"$run_ref\""/>" >> "${alias}.xml"
  echo "        <ANALYSIS_TYPE>" >> "${alias}.xml"
  echo "            <PROCESSED_READS/>" >> "${alias}.xml"
  echo "        </ANALYSIS_TYPE>" >> "${alias}.xml"
  echo "        <FILES>" >> "${alias}.xml"
  echo "            <FILE filename="\"$file_names\"" filetype="\"bam\""" >> "${alias}.xml"
  echo "                checksum_method="\"$checksum_method\"" checksum="\"$checksum\""/>" >> "${alias}.xml"
  echo "        </FILES>" >> "${alias}.xml"
  echo "        <ANALYSIS_ATTRIBUTES>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Project</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$project</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Assay Type</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$assay_type</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Analysis protocol</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$analysis_protocol</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Analysis code</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$analysis_code</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Reference genome</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$reference_genome</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Analysis center</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$analysis_center</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "            <ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "                <TAG>Analysis date</TAG>" >> "${alias}.xml"
  echo "                <VALUE>$analysis_date</VALUE>" >> "${alias}.xml"
  echo "            </ANALYSIS_ATTRIBUTE>" >> "${alias}.xml"
  echo "        </ANALYSIS_ATTRIBUTES>" >> "${alias}.xml"
  echo "    </ANALYSIS>" >> "${alias}.xml"
  echo "</ANALYSIS_SET>" >> "${alias}.xml"

done < "$metadata"


# Create analysis submission XML (only once)

echo "<SUBMISSION>" >> submission.xml
echo "   <ACTIONS>" >> submission.xml
echo "      <ACTION>" >> submission.xml
echo "         <ADD/>" >> submission.xml
echo "      </ACTION>" >> submission.xml
echo "   </ACTIONS>" >> submission.xml
echo "</SUBMISSION>" >> submission.xml