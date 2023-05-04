### Script for grabbing input for the 1_preprocessing.snakefile from the EHI AirTable.
### This pulls EHI numbers from the 'PR_preprocessing' table.
### Raphael Eisenhofer 4/2023

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--prb', required=True, help='preprocessing batch number - e.g. "PRB0001"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblJfLRU2FIVz37Y1'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{PR_batch}} = '{args.prb}'",
    'fields': ['EHI_number']

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'prb_input.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')

    for record in records:
        # Get the values of the EHI_number lookup fields
        # EHI_number are lookups in the AB_assembly_binning table, so we perform further API requests.
        ehi_number_id = record['fields']['EHI_number'][0]

        # Make requests to retrieve the linked records
        ehi_number_response = requests.get(f"{AIRTABLE_API_ENDPOINT}/{ehi_number_id}", headers=headers)

        # Extract the values of the linked fields from the linked records
        ehi_number_value = ehi_number_response.json()['fields']['EHI_number']

        # Write the row to the TSV file
        row = [ehi_number_value]
        writer.writerow(row)