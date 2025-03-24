### Script for grabbing input for the EHIO cataloging pipeline
### This pulls preprocessed read URLs from the 'PR_preprocessing' table.
### Raphael Eisenhofer 3/2025

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--asb', required=True, help='assembly batch number - e.g. "ABB0001"')
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
    'filterByFormula': f"{{Assembly batches}} = '{args.asb}'",

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
    writer.writerow(['EHI_number', 'metagenomic_bases', 'r1', 'r2'])

    offset = None  # Initialize offset
    while True:
        # Include offset in query parameters if it exists
        if offset:
            query_params['offset'] = offset

        # Make the request
        response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)
        if response.status_code != 200:
            raise Exception(f"Error: {response.status_code}, {response.text}")

        # Parse response JSON
        data = response.json()

        # Process records
        for record in data.get('records', []):
            fields = record['fields']
            ehi = fields.get('EHI_plaintext')
            metagenomic_bases = fields.get('metagenomic_bases')
            url1 = fields.get('URL_meta1', '')
            url2 = fields.get('URL_meta2', '')

            writer.writerow([ehi, metagenomic_bases, url1, url2])

        if 'offset' in data:
            offset = data['offset']
        else:
            break