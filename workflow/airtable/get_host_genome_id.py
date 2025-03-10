### Get EHI host genome ID
## Raphael Eisenhofer 3/2025

import argparse
import requests
import csv
import json

#Add input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--prb', required=True, help='PRB ID')
args = parser.parse_args()
prb_value = args.prb

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tbl1l5mAkF9nVtisn'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

    #Get record ID from AirTable
params = {
        'filterByFormula': f"Code = '{args.prb}'",
        'maxRecords': 1
    }
response = requests.get(url, headers=headers, params=params)

record_id = response.json().get('records')[0].get('id')

# Set up the output TSV file
output_file_path = 'host_genome.tsv'

# Check if the request was successful
if response.status_code == 200:
    records = response.json().get('records')
    if records:
        record = records[0]
        record_id = record.get('id')
        
        # Extract the value from airtable
        host_genome_value = record.get('fields', {}).get('reference_genome_plain')
        
        if host_genome_value:
            # Set up the output TSV file
            output_file_path = 'host_genome.tsv'
            
            with open(output_file_path, 'w', newline='') as tsvfile:
                writer = csv.writer(tsvfile, delimiter='\t')
                
                # Write the header and the value to the TSV file
                writer.writerow([host_genome_value])
                
            print(f"Data written to {output_file_path}")
        else:
            print("No value found in 'Reference genome'")
    else:
        print("No records found matching the code.")
else:
    print(f"Failed to fetch data from Airtable. Status code: {response.status_code}")