## Raphael Eisenhofer 2023

import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--mags', required=True, help='Path to list of dereplicated MAGs')
parser.add_argument('--dmb', required=True, help='DMB number')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblvPsCiNanM7NeQp'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}


# Read the entire content of the --mags file
with open(args.mags, 'r') as file:
    mags_content = file.read()

mags_content = mags_content.replace('\n', ',')

# Get the record ID for the row based on the value in the 'DM_batch'
params = {
    'filterByFormula': f"{{Code}} = '{args.dmb}'",
    'maxRecords': 1
}
response = requests.get(url, headers=headers, params=params)
data = response.json()
record_id = data['records'][0]['id']

# Set the cell data you want to update
data = {
    'fields': {
        'dereplicated_mags': mags_content,
    }
}

# Send a PATCH request to update the record
response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

# Print the response status code
print(response.status_code)