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
    'filterByFormula': f"FIND('{args.dmb}', {{DM_batch}})",

}

initial_offset = 1
all_records = []

# Loop until there are no more pages
while initial_offset != 0:
  print(initial_offset)
  response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)
  data = response.json()
  records = data['records']    
  all_records.extend(records)  
  if 'offset' in data:
    initial_offset = data['offset']
    query_params['offset'] = initial_offset
  else:
    initial_offset = 0
    break

# Set up the output TSV file
output_file_path = 'read_input.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['PR_batch', 'EHI_number'])

    offset = None
    while True:
        # Update the query parameters with the offset if it exists
        if offset:
            query_params['offset'] = offset

        for record in all_records:
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

        if 'offset' in data:
            offset = data['offset']
        else:
            break