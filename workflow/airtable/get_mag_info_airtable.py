## Raphael Eisenhofer 2023

import argparse
import csv
import json
import requests

# Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--dmb', required=True, help='DMB number')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblMzd3oyaJhdeQcs'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"AND(FIND('{args.dmb}', {{DM_batch}}), {{dereplicated}} = 'true')",
    'pageSize': 100  # Set the page size to 100 records per request
}

# Set up the output TSV file
output_file_path = f'{args.dmb}_mag_info.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['genome', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'completeness', 'contamination', 'mag_size'])

    offset = None
    while True:
        # Update the query parameters with the offset if it exists
        if offset:
            query_params['offset'] = offset

        # Make the request to get the records
        response = requests.get(url, params=query_params, headers=headers)

        if response.status_code != 200:
            print('Error: Failed to retrieve records from the Airtable API.')
            print('Status Code:', response.status_code)
            break

        # Convert the response to a JSON object
        data = response.json()

        # Extract the records from the JSON object
        records = data['records']

        # Loop through the records and write each one to the CSV file
        for record in records:
            record_id = record['id']
            record_response = requests.get(f"{url}/{record_id}", headers=headers)

            # Extract the values of the completeness and contamination fields
            mag_name = record_response.json()['fields'].get('mag_name', '')
            domain = record_response.json()['fields'].get('domain', '')
            phylum = record_response.json()['fields'].get('phylum', '')
            phyclass = record_response.json()['fields'].get('class', '')
            order = record_response.json()['fields'].get('order', '')
            family = record_response.json()['fields'].get('family', '')
            genus = record_response.json()['fields'].get('genus', '')
            species = record_response.json()['fields'].get('species', '')
            completeness = record_response.json()['fields'].get('completeness', '')
            contamination = record_response.json()['fields'].get('contamination', '')
            mag_size = record_response.json()['fields'].get('size', '')

            # Write the row to the CSV file
            row = [mag_name, domain, phylum, phyclass, order, family, genus, species, completeness, contamination, mag_size]
            writer.writerow(row)

        if 'offset' in data:
            offset = data['offset']
        else:
            break
