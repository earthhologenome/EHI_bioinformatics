### Script for experiment checklist for ENA registration
### This pulls sample numbers from the 'Samples' table.
### Raphael Eisenhofer 6/2024

import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--ehi', required=True, help='EHI ID - e.g. "EHI01280"')
#parser.add_argument('--sample', required=True, help='Sample code - e.g. "ARS50"')
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
    'filterByFormula': f"{{EHI_number}} = '{args.ehi}'",

}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']


column_names = [
    "alias", "title", "study_alias", "sample_alias", "design_description", "library_name", "insert_size", "library_layout", "library_strategy",
    "library_source", "library_selection", "platform", "instrument_model", "forward_file_name", "reverse_file_name" 
]


# Set up the output TSV file
output_file_path = f'{args.ehi}_experiment_checklist.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(column_names)

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
            alias_value = record_response.json()['fields'].get('sample_alias', '')
            title_value = record_response.json()['fields'].get('EHI_number', '')
            study_alias_value = record_response.json()['fields'].get('study_alias', '')
            sample_alias_value = record_response.json()['fields'].get('sample_alias', '')
            design_description_value = record_response.json()['fields'].get('design_description', '')
            library_name_value = record_response.json()['fields'].get('library_name', '')
            insert_size_value = record_response.json()['fields'].get('insert_size', '')
            library_layout_value = record_response.json()['fields'].get('library_layout', '')
            library_strategy_value = record_response.json()['fields'].get('library_strategy', '')
            library_source_value = record_response.json()['fields'].get('library_source', '')
            library_selection_value = record_response.json()['fields'].get('library_selection', '')
            platform_value = record_response.json()['fields'].get('platform', '')
            instrument_model_value = record_response.json()['fields'].get('instrument_model', '')
            forward_file_name_value = f'{args.ehi}_raw_1.fq.gz'
            reverse_file_name_value = f'{args.ehi}_raw_2.fq.gz'

            # Write the row to the TSV file
            row = [
                alias_value, title_value, study_alias_value, sample_alias_value, design_description_value, library_name_value, insert_size_value, library_layout_value,
                library_strategy_value, library_source_value, library_selection_value, platform_value, instrument_model_value,
                forward_file_name_value, reverse_file_name_value
            ]
            writer.writerow(row)

        if 'offset' in data:
            offset = data['offset']
        else:
            break