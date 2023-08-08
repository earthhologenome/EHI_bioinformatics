import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--report', required=True, help='Path to report file')
parser.add_argument('--prb', required=True, help='PR batch number')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblJfLRU2FIVz37Y1'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

#Read in the TSV file using pandas
df = pd.read_csv(args.report, sep='\t')

# Loop through each row in the dataframe
for i, row in df.iterrows():
    # Get the record ID for the row based on the value in the 'ehi_number' column and 'pr_batch'
    params = {
        'filterByFormula': f"AND({{EHI_number}} = '{row['EHI_number']}', {{PR_Batch}} = '{args.prb}')",        
        'maxRecords': 1
    }
    response = requests.get(url, headers=headers, params=params)
    data = response.json()
    record_id = data['records'][0]['id']

    # Set the cell data you want to update
    data = {
        'fields': {
            'reads_pre_fastp': int(row['reads_pre_fastp']),
            'reads_post_fastp': int(row['reads_post_fastp']),
            'bases_pre_fastp': int(row['bases_pre_fastp']),
            'bases_post_fastp': int(row['bases_post_fastp']),
            'adapter_trimmed_reads': int(row['adapter_trimmed_reads']),
            'adapter_trimmed_bases': int(row['adapter_trimmed_bases']),
            'host_reads': int(row['host_reads']),
            'bacterial_archaeal_bases': int(row['bacterial_archaeal_bases']),
            'metagenomic_bases': int(row['metagenomic_bases']),
            'singlem_fraction': float(row['singlem_fraction'].strip('%'))/100,
            'kappa': float(row['kappa']),
            'C': float(row['C']),
            'LR': int(row['LR']),
            'modelR': float(row['modelR']),
            'LRstar': int(row['LRstar']),
            'diversity': float(row['diversity']),
            'host_duplicate_fraction': float(row['host_duplicate_fraction'])
        }
    }

    # Send a PATCH request to update the record
    response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)