import argparse
import csv
import pandas as pd
import requests
import json

# Load CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--ehs', required=True, help='EHS number')
parser.add_argument('--ehi', required=True, help='EHI number')
args = parser.parse_args()

# Read the API key from the config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
base_url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblMzd3oyaJhdeQcs'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"AND(FIND('{args.ehs}', {{Submission}}), NOT({{phylum}} = 'nan'))",
    'pageSize': 100  # Set the page size to 100 records per request
}

# Set up the output TSV file
output_file_path = f'{args.ehi}_mag_checklist_temp.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(['mag_id', 'mag_name', 'metagenomic source', 'assembly software', 
    'binning software', 'binning parameters', 'assembly quality', 'completeness software',
    'completeness score', 'contamination score', 'taxonomic identity marker', 'taxonomic classification'])

    offset = None
    while True:
        # Update the query parameters with the offset if it exists
        if offset:
            query_params['offset'] = offset

        # Make the request to get the records
        response = requests.get(base_url, params=query_params, headers=headers)

        # Convert the response to a JSON object
        data = response.json()

        # Extract the records from the JSON object
        records = data['records']

        # Loop through the records and write each one to the CSV file
        for record in records:
            # Get the ID of the record
            record_id = record['id']

            # Make a request to get the record data
            record_response = requests.get(f"{base_url}/{record_id}", headers=headers)

            # Extract the values of the completeness and contamination fields
            mag_id = record_response.json()['fields'].get('ID', '')
            mag_name = record_response.json()['fields'].get('mag_name', '')
            metagenomic_source = record_response.json()['fields'].get('metagenomic source', '')
            assembly_software = record_response.json()['fields'].get('assembly software', '')
            binning_software = record_response.json()['fields'].get('binning software', '')
            binning_parameters = record_response.json()['fields'].get('binning parameters', '')
            assembly_quality = record_response.json()['fields'].get('assembly quality', '')
            completeness_software = record_response.json()['fields'].get('completeness software', '')
            completeness = record_response.json()['fields'].get('completeness', '')
            contamination = record_response.json()['fields'].get('contamination', '')
            taxonomic_identity_marker = record_response.json()['fields'].get('taxonomic identity marker', '')
            taxonomic_classification = record_response.json()['fields'].get('taxonomic classification', '')

            # Write the row to the CSV file
            row = [mag_id, mag_name, metagenomic_source, assembly_software, binning_software,
            binning_parameters, assembly_quality, completeness_software, completeness, contamination,
            taxonomic_identity_marker, taxonomic_classification]
            writer.writerow(row)

        # Check if there are more records to retrieve
        if 'offset' in data:
            offset = data['offset']
        else:
            break



output_file_path = f'{args.ehi}_mag_urls.tsv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')

    offset = None
    while True:
        # Update the query parameters with the offset if it exists
        if offset:
            query_params['offset'] = offset

        # Make the request to get the records
        response = requests.get(base_url, params=query_params, headers=headers)

        # Convert the response to a JSON object
        data = response.json()

        # Extract the records from the JSON object
        records = data['records']

        # Loop through the records and write each one to the CSV file
        for record in records:
            # Get the ID of the record
            record_id = record['id']

            # Make a request to get the record data
            record_response = requests.get(f"{base_url}/{record_id}", headers=headers)

            # Extract the values of the completeness and contamination fields
            MAG_url = record_response.json()['fields'].get('MAG_url', '')

            # Write the row to the CSV file
            row = [MAG_url]
            writer.writerow(row)

        # Check if there are more records to retrieve
        if 'offset' in data:
            offset = data['offset']
        else:
            break