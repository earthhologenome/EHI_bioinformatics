### Automatically update AirTable from Mjolnir (with EHS completion status)
## Raphael Eisenhofer 6/2024

import argparse
import requests
import json

#Add input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ehs', required=True, help='EHS number')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblDxjHlzqMfwxvPN'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

    #Get record ID from AirTable
params = {
        'filterByFormula': f"ID = '{args.ehs}'",
        'maxRecords': 1
    }
response = requests.get(url, headers=headers, params=params)
record_id = response.json().get('records')[0].get('id')


    #Change the value in the AirTable
data = {
        'fields': {
            'Status': 'Running',
        }
    }

# Send a PATCH request to update the record in the AirTable
response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

# Did it work?
print(response.status_code)