### Script for grabbing input for the 5_ena_upload.snakefile from the EHI AirTable.
### This pulls EHI numbers from the 'SE (Samples)' table.
### Raphael Eisenhofer 6/2024

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--ehs', required=True, help='preprocessing batch number - e.g. "EHS001"')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblaMWLkBUn2g5gRR'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{Submissions}} = '{args.ehs}'"

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'ehi_numbers.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['EHI', 'SAMPLE', 'r1', 'r2'])

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
            ehi_number_value = fields.get('EHI_number', '')
            sample_number_value = fields.get('sample_alias', '')
            forward_url_value = fields.get('forward_url', '')
            reverse_url_value = fields.get('reverse_url', '')

            # Write to file
            writer.writerow([ehi_number_value, sample_number_value, forward_url_value, reverse_url_value])

        # Check for next page
        offset = data.get('offset')
        if not offset:
            break
