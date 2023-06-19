## Raphael Eisenhofer 2023

import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--report', required=True, help='Path to report file')
parser.add_argument('--asb', required=True, help='ASB number')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblG6ZIvkYN844I97'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

#Read in the TSV file using pandas
df = pd.read_csv(args.report, sep='\t')

# Loop through each row in the dataframe
for i, row in df.iterrows():
    # Get the record ID for the row based on the value in the 'EHA_number' column and 'AB_batch'
    params = {
        'filterByFormula': f"AND({{Assembly_code}} = '{row['EHA_number']}', {{EHI_number_api}} = '{row['EHI_number']}')",        
        'maxRecords': 1
    }
    response = requests.get(url, headers=headers, params=params)
    data = response.json()
    record_id = data['records'][0]['id']

    # Set the cell data you want to update
    data = {
        'fields': {
            'assembly_length': int(row['assembly_length']),
            'N50': int(row['N50']),
            'L50': int(row['L50']),
            'num_contigs': int(row['num_contigs']),
            'largest_contig': int(row['largest_contig']),
            'num_bins': int(row['num_bins']),
            'assembly_mapping_percent': int(row['assembly_mapping_percent'])
        }
    }

    # Send a PATCH request to update the record
    response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)