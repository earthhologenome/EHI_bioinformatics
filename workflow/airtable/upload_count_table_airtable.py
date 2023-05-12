https://chinarajames.com/how-to-upload-images-to-airtable-using-their-api/

https://sid.erda.dk/share_redirect/BaMZodj9sA/


import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--table', required=True, help='Path to count table')
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

#Read file filter by DMB
params = {
    'filterByFormula': f'{{Code}} = "{args.dmb}"'
}

# Send a GET request to Airtable to retrieve the record ID
response = requests.get(url, headers=headers, params=params)
data = response.json()
if len(data['records']) == 0:
    print(f'No record found with Code "{args.dmb}"')
    exit()
record_id = data['records'][0]['id']

# Open the gzip file and read the contents as bytes
with open(args.table, 'rb') as f:
    contents = f.read()

# Set the data for the request
data = {
    'fields': {
        'count_table': [
            {
                'url': f'https://api.airtable.com/v0/{url}/{record_id}/attachmentName',
                'filename': args.table,
                'size': len(contents),
                'type': 'application/gzip'
            }
        ]
    }
}

# Set the URL for the PATCH request to update the record
patch_url = f'https://api.airtable.com/v0/{url}/{record_id}'

# Send the PATCH request to Airtable to update the record with the attachment
response = requests.patch(patch_url, headers=headers, files={'AttachmentColumn': (args.table, contents)})