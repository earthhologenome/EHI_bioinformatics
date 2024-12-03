### API call for getting ENA sample accession
### This pulls sample numbers from the 'Samples' table.
### Raphael Eisenhofer 12/2024

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--ehi', required=True, help='EHI ID - e.g. "EHI01280"')
parser.add_argument('--sample', required=True, help='Sample code - e.g. "ARS50"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblW03Z3DcjRdEkoS'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{Code}} = '{args.sample}'",

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = f'{args.ehi}_ENA_sample_accession.txt'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')

    offset = None
    while True:
        # Update the query parameters with the offset if it exists
        if offset:
            query_params['offset'] = offset

        for record in records:
            # Get the values of the PR_batch and EHI_number lookup fields
            record_id = record['id']

            # Make requests to retrieve the linked records
            record_response = requests.get(f"{AIRTABLE_API_ENDPOINT}/{record_id}", headers=headers)

            # Extract the values of the linked fields from the linked records
            ena_sample_accession_value = record_response.json()['fields'].get('ENA_sample_accession', '')

            # Write the row to the TSV file
            row = [
                ena_sample_accession_value
            ]
            writer.writerow(row)

        if 'offset' in data:
            offset = data['offset']
        else:
            break