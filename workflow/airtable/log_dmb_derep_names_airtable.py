# Raphael Eisenhofer 7/2023

import argparse
import pandas as pd
import requests
import json

# Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--mags', required=True, help='Path to list of dereplicated MAGs')
parser.add_argument('--dmb', required=True, help='DMB code')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
dmburl = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblvPsCiNanM7NeQp'
magurl = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblMzd3oyaJhdeQcs'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Search for the args.dmb value in the dmburl table to get the corresponding record ID
params = {
    'filterByFormula': f"{{Code}} = '{args.dmb}'",
}
response = requests.get(dmburl, headers=headers, params=params)
data = response.json()
dmb_record_id = data['records'][0]['id']


# Read the contents of the TSV file
with open(args.mags, 'r') as file:
    mags_content = file.read().splitlines()

# Iterate over the MAGs and update 'dereplicated_test' field for matching records
for mag_name in mags_content:
    params = {
        'maxRecords': 1,
        'filterByFormula': f"{{mag_name}} = '{mag_name}'"
    }
    response = requests.get(magurl, headers=headers, params=params)
    data = response.json()
    if data['records']:
        record_id = data['records'][0]['id']

        # Retrieve the existing values of 'dereplicated_test' field
        existing_values = data['records'][0]['fields'].get('dereplicated', [])

        # Append the new value to the existing values
        existing_values.append(dmb_record_id)

        data = {
            'fields': {
                'dereplicated': existing_values  # Use the retrieved dmb_record_id here
            }
        }
        response = requests.patch(f'{magurl}/{record_id}', headers=headers, data=json.dumps(data))
        print(f"Updated record {record_id}. Response status code: {response.status_code}")
    else:
        print(f"MAG '{mag_name}' not found in the 'magurl' table.")
