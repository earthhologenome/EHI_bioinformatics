## Raphael Eisenhofer 2023

import argparse
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--table', required=True, help='URL to count table')
parser.add_argument('--tree', required=True, help='URL to the tree')
parser.add_argument('--taxonomy', required=True, help='URL to the taxonomy')
parser.add_argument('--dmb', required=True, help='DMB identifier')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# set up the API endpoint and headers
api_endpoint = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblvPsCiNanM7NeQp'
headers = {'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'}

# get the records from the Airtable where 'Code' equals the value of --dmb
params = {'filterByFormula': f"{{Code}}='{args.dmb}'"}
response = requests.get(api_endpoint, headers=headers, params=params).json()

record_id = response['records'][0]['id']  # get the record ID
update_data = {
    "fields": {
        "phylogenetic_tree": args.tree,
        "taxonomy_table": args.taxonomy,
        "count_table": args.table
    }
}
update_url = f"{api_endpoint}/{record_id}"
update_response = requests.patch(update_url, headers=headers, json=update_data)

if update_response.status_code == 200:
    print("Record updated successfully.")
else:
    print(f"Failed to update record. Error code: {update_response.status_code}")
