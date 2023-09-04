import argparse
import csv
import pandas as pd
import requests
import json

# Load CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--dmb', required=True, help='DMB number')
args = parser.parse_args()

# Read the API key from the config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
base_url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblvPsCiNanM7NeQp'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"FIND('{args.dmb}', {{Code}})",
    'pageSize': 100  # Set the page size to 100 records per request
}

# Set up the output TSV file
output_file_path = 'num_mags.csv'


response = requests.get(base_url, headers=headers, params=query_params)

# Check if the request was successful (status code 200)
if response.status_code == 200:
    data = response.json()
    
    # Extract the 'number_of_input_MAGs' value from the first record (assuming you only need one)
    if data.get('records'):
        num_mags = data['records'][0]['fields']['number_of_input_MAGs']
        
        # Write the numerical value to 'num_mags.csv'
        with open('num_mags.csv', 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(['number_of_input_MAGs'])  # Write header
            csv_writer.writerow([num_mags])  # Write the numerical value
    
else:
    print(f"Failed to retrieve data. Status code: {response.status_code}")
