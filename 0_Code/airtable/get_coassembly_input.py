### Script for grabbing input for the 2_Assembly_mag.snakefile from the EHI AirTable.
### This pulls EHA, PRB, and EHI numbers from the 'AB_assembly_binning' table.
### Raphael Eisenhofer 4/2023

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--abb', required=True, help='Assembly batch number - e.g. "ABB0001"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblG6ZIvkYN844I97'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{AB_batch}} = '{args.abb}'",
    'fields': ['ID', 'PR_batch', 'EHI_number', 'Assembly_code']

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'asb_input.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['PR_batch', 'EHI_number', 'Assembly_code'])

    for record in records:
        # Get the values of the PR_batch and EHI_number lookup fields
        # PR_batch and EHI_number are lookups in the AB_assembly_binning table, so we perform further API requests.
        pr_batch_id = record['fields']['PR_batch'][0]
        ehi_number_id = record['fields']['EHI_number'][0]

        # Make requests to retrieve the linked records
        pr_batch_response = requests.get(f"{AIRTABLE_API_ENDPOINT}/{pr_batch_id}", headers=headers)
        ehi_number_response = requests.get(f"{AIRTABLE_API_ENDPOINT}/{ehi_number_id}", headers=headers)

        # Extract the values of the linked fields from the linked records
        # Note that PR_batch is called 'Code' in the PR_batch table.
        pr_batch_value = pr_batch_response.json()['fields']['Code']
        ehi_number_value = ehi_number_response.json()['fields']['EHI_number']

        # Write the row to the TSV file
        row = [pr_batch_value, ehi_number_value, record['fields']['Assembly_code'][0]]
        writer.writerow(row)