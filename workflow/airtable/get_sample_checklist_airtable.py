### Script for sample checklist for ENA registration
### This pulls sample numbers from the 'Samples' table.
### Raphael Eisenhofer 6/2024

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


column_names = [
    "alias", "sample_alias", "title", "sample_title", "sample_description", "host subject id", "host common name", "host taxid", "host lifestage",
    "project name", "taxon_id", "scientific_name", "collection date",
    "geographic location (country and/or sea)", "geographic location (latitude)", "geographic location (longitude)", "broad-scale environmental context", "local environmental context",
    "environmental medium"
]


# Set up the output TSV file
output_file_path = f'{args.ehi}_sample_checklist.tsv'

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
            alias_value = "TEMP"
            sample_alias_value = "TEMP"
            title_value = "TEMP"
            sample_title_value = "TEMP"
            sample_description_value = record_response.json()['fields'].get('sample_description', '')
            host_subject_id_value = record_response.json()['fields'].get('host subject id', '')
            host_common_name_value = record_response.json()['fields'].get('host common name', '')
            host_taxid_value = record_response.json()['fields'].get('host taxid', '')
            host_lifestage_value = record_response.json()['fields'].get('host lifestage', '')
            project_name_value = record_response.json()['fields'].get('project name', '')
            taxon_id_value = record_response.json()['fields'].get('taxon_id', '')
            scientific_name_value = record_response.json()['fields'].get('scientific_name', '')
            collection_date_value = record_response.json()['fields'].get('collection date', '')
            geographic_location_value = record_response.json()['fields'].get('geographic location (country and/or sea)', '')
            geo_lat_value = record_response.json()['fields'].get('geographic location (latitude)', '')
            geo_long_value = record_response.json()['fields'].get('geographic location (longitude)', '')
            broad_env_value = record_response.json()['fields'].get('broad-scale environmental context', '')
            local_env_value = record_response.json()['fields'].get('local environmental context', '')
            env_medium_value = record_response.json()['fields'].get('environmental medium', '')

            # Write the row to the TSV file
            row = [
                alias_value, sample_alias_value, title_value, sample_title_value, sample_description_value, host_subject_id_value, host_common_name_value, host_taxid_value, host_lifestage_value,
                project_name_value, taxon_id_value, scientific_name_value, collection_date_value,
                geographic_location_value, geo_lat_value, geo_long_value, broad_env_value, local_env_value,
                env_medium_value
            ]
            writer.writerow(row)

        if 'offset' in data:
            offset = data['offset']
        else:
            break