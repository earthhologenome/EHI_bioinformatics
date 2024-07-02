import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--ehi', required=True, help='EHI number')
parser.add_argument('--exp_acc', required=True, help='ENA experiment accession')
parser.add_argument('--run_acc', required=True, help='ENA run accession')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblaMWLkBUn2g5gRR'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

params = {
    'filterByFormula': f"{{EHI_number}} = '{args.ehi}'",

    }
response = requests.get(url, headers=headers, params=params)
data = response.json()
record_id = data['records'][0]['id']

    # Set the cell data you want to update
data = {
        'fields': {
            'ENA_experiment_accession': {args.exp_acc},
            'ENA_run_accession': {args.run_acc},
        }
    }

    # Send a PATCH request to update the record
response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

    # Print the response status code
print(response.status_code)