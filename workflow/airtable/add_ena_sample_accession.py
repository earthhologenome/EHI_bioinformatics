import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--sample', required=True, help='EHI sample code')
parser.add_argument('--sample_acc', required=True, help='ENA sample accession')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblW03Z3DcjRdEkoS'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}


    # Get the record ID for the row based on the value in the 'ehi_number' column and 'pr_batch'
params = {
    'filterByFormula': f"{{Code}} = '{args.sample}'",

    }
response = requests.get(url, headers=headers, params=params)
data = response.json()
record_id = data['records'][0]['id']

    # Set the cell data you want to update
data = {
        'fields': {
            'ENA_sample_accession': {args.sample_acc},
        }
    }

    # Send a PATCH request to update the record
response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

    # Print the response status code
print(response.status_code)