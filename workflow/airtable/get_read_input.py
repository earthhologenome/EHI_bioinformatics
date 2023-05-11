### Script for grabbing reads for the 4_dereplication_mapping.snakefile from the EHI AirTable.
### This pulls EHI numbers from the 'PR_preprocessing' table.
### Raphael Eisenhofer 5/2023

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--dmb', required=True, help='DMB batch number - e.g. "DMB0001"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblkmD4y0fZ6P8XRf'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{DM_batch}} = '{args.dmb}'",

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'read_input.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['PR_batch', 'EHI_number'])

    for record in records:
        # Get the values of the PR_batch and EHI_number lookup fields
        record_id = record['id']

        # Make requests to retrieve the linked records
        record_response = requests.get(f"{AIRTABLE_API_ENDPOINT}/{record_id}", headers=headers)

        # Extract the values of the linked fields from the linked records
        pr_batch_value = record_response.json()['fields'].get('PR_Batch', '')
        ehi_number_value = record_response.json()['fields'].get('EHI_number', '')

        # Write the row to the TSV file
        row = [pr_batch_value, ehi_number_value]
        writer.writerow(row)