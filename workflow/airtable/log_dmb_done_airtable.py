### Automatically update AirTable from Mjolnir (with PRB completion status)
## Raphael Eisenhofer 3/2023

import argparse
import requests
import json

#Add input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--code', required=True, help='DMB ID')
args = parser.parse_args()
code_value = args.code

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

    #Get record ID from AirTable
params = {
        'filterByFormula': f"Code = '{args.code}'",
        'maxRecords': 1
    }
response = requests.get(url, headers=headers, params=params)
record_id = response.json().get('records')[0].get('id')


    #Change the value in the AirTable
data = {
        'fields': {
            'Status': 'Done',
        }
    }

# Send a PATCH request to update the record in the AirTable
response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

# Did it work?
print(response.status_code)