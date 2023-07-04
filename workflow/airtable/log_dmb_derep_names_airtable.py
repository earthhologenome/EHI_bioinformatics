## Raphael Eisenhofer 2023

import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--mags', required=True, help='Path to list of dereplicated MAGs')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblMzd3oyaJhdeQcs'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}


# Read the contents of the TSV file
with open(args.mags, 'r') as file:
    mags_content = file.read().splitlines()

# Iterate over the MAGs and update 'dereplicated' column for matching records
for mag_name in mags_content:
    params = {
        'maxRecords': 1,
        'filterByFormula': f"{{mag_name}} = '{mag_name}'"
    }
    response = requests.get(url, headers=headers, params=params)
    data = response.json()
    if data['records']:
        record_id = data['records'][0]['id']
        data = {
            'fields': {
                'dereplicated': 'true',
            }
        }
        response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))
        print(f"Updated record {record_id}. Response status code: {response.status_code}")