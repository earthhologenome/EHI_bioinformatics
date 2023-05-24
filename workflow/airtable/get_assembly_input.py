### Script for grabbing input for the 2_Assembly_mag.snakefile from the EHI AirTable.
### This pulls EHA, PRB, and EHI numbers from the 'AB_assembly_binning' table.
### Raphael Eisenhofer 4/2023

import argparse
import json
import requests
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--abb', required=True, help='Assembly batch number - e.g. "ABB0001"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_ASB = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblG6ZIvkYN844I97'
AIRTABLE_PPR = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblJfLRU2FIVz37Y1'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters for assembly table to filter the records
query_params = {
    'filterByFormula': f"{{AB_batch}} = '{args.abb}'",
    'fields': ['ID', 'PR_batch', 'EHI_number', 'Assembly_code']
}

# Make the request to get the records from assembly table
response = requests.get(AIRTABLE_ASB, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object for assembly table
asb_records = data['records']

# Set up the output TSV file
output_file_path = 'asb_input.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['PR_batch', 'EHI_number', 'Assembly_code', 'metagenomic_bases', 'singlem_fraction', 'diversity', 'nonpareil_estimated_coverage'])

    # Process records from assembly table
    for record in asb_records:
        # Get the values of the PR_batch and EHI_number lookup fields
        pr_batch_id = record['fields']['PR_batch'][0]
        ehi_number_id = record['fields']['EHI_number'][0]

        # Make requests to retrieve the linked records
        pr_batch_response = requests.get(f"{AIRTABLE_ASB}/{pr_batch_id}", headers=headers)
        ehi_number_response = requests.get(f"{AIRTABLE_ASB}/{ehi_number_id}", headers=headers)

        # Extract the values of the linked fields from the linked records
        pr_batch_value = pr_batch_response.json()['fields']['Code']
        ehi_number_value = ehi_number_response.json()['fields']['EHI_number']

        # Set up the query parameters for Table 2 based on the 'EHI_number' value
        ppr_query_params = {
            'filterByFormula': f"{{EHI_number}} = '{ehi_number_value}'",
            'fields': ['metagenomic_bases', 'singlem_fraction', 'diversity', 'nonpareil_estimated_coverage']
        }

        # Make the request to get the records from Table 2
        ppr_response = requests.get(AIRTABLE_PPR, params=ppr_query_params, headers=headers)

        # Convert the response to a JSON object
        ppr_data = ppr_response.json()

        # Extract the records from the JSON object for Table 2
        ppr_records = ppr_data['records']

        # Process records from Table 2
        for ppr_record in ppr_records:
            # Extract the values from Table 2 and perform further processing
            metagenomic_bases_value = ppr_record['fields']['metagenomic_bases']
            singlem_fraction_value = ppr_record['fields']['singlem_fraction']
            diversity_value = ppr_record['fields']['diversity']
            nonpareil_estimated_coverage_value = ppr_record['fields']['nonpareil_estimated_coverage'].replace('%', '')

            # Write the row to the TSV file
            row = [pr_batch_value, ehi_number_value, record['fields']['Assembly_code'], metagenomic_bases_value, singlem_fraction_value, diversity_value, nonpareil_estimated_coverage_value]
            writer.writerow(row)
