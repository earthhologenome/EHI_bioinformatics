## Raphael Eisenhofer 2023

import argparse
import pandas as pd
import requests
import json

# Load CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--report', required=True, help='Path to genome report file')
args = parser.parse_args()

# Read the API key from the config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
base_url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV'
table_name = 'tblMzd3oyaJhdeQcs'
url = f'{base_url}/{table_name}'

headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Read the TSV file using pandas
df = pd.read_csv(args.report, sep='\t')

# Loop through each row in the dataframe
for i, row in df.iterrows():
    # Check if the 'mag_name' already exists in Airtable
    params = {
        'filterByFormula': f"{{mag_name}} = '{row['mag_name']}'"
    }
    existing_records_response = requests.get(url, headers=headers, params=params)
    existing_records = existing_records_response.json().get('records', [])
    existing_record = next(iter(existing_records), None)

    # Set the cell data you want to update
    data = {
        'fields': {
            'mag_name': str(row['mag_name']),
            'eha_number': str(row['eha_number']),
            'GTDB_version': str(row['GTDB_version']),
            'domain': str(row['domain']),
            'phylum': str(row['phylum']),
            'class': str(row['class']),
            'order': str(row['order']),
            'family': str(row['family']),
            'genus': str(row['genus']),
            'species': str(row['species']),
            'fastani_ani': float(row['fastani_ani']),
            'closest_placement_af': float(row['closest_placement_af']),
            'closest_placement_ani': float(row['closest_placement_ani']),
            'completeness': float(row['completeness']),
            'contamination': float(row['contamination']),
            'size': int(row['size']),
            'GC': float(row['GC']),
            'coding_density': float(row['coding_density']),
            'N50': int(row['N50']),
            'contigs': int(row['contigs'])
        }
    }

    if existing_record:
        # Update existing record
        record_id = existing_record['id']
        update_url = f'{url}/{record_id}'
        response = requests.patch(update_url, headers=headers, data=json.dumps(data))
    else:
        # Create a new record
        response = requests.post(url, headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)
