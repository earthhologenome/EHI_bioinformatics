
import argparse
import pandas as pd
import requests
import json

# Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--samples', required=True, help='Path to samples file')
parser.add_argument('--dmb', required=True, help='DMB number')
args = parser.parse_args()

# Read the API key from the config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set variables
url = 'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tblJfLRU2FIVz37Y1'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

# Read in the TSV file using pandas
df = pd.read_csv(args.samples, sep='\t')

# Create an empty list to store the filtered records
filtered_records = []

# Loop through each row in the dataframe
for i, row in df.iterrows():
    # Get the record ID for the row based on the value in the 'EHI_number' column and 'PR_batch'
    params = {
        'filterByFormula': f"AND({{EHI_number}} = '{row['EHI_number']}', {{PR_Batch}} = '{row['PR_batch']}')",
        'maxRecords': 1
    }
    response = requests.get(url, headers=headers, params=params)
    data = response.json()

    # Check if any records are found
    if 'records' in data and len(data['records']) > 0:
        record_fields = data['records'][0]['fields']
        filtered_records.append(record_fields)



# Specify the desired columns to be included in the new TSV file
desired_columns = ['EHI_plaintext', 'sample_code', 'species', 'region', 'sample_type', 'order', 'sex', 'country', 'latitude', 'longitude', 'singlem_fraction', 'metagenomic_bases', 'host_bases', 'bases_lost_fastp_percent', 'diversity', 'C']

# Create a new DataFrame from the filtered records
filtered_records_df = pd.DataFrame(filtered_records)

# Select only the desired columns
filtered_records_df = filtered_records_df[desired_columns]

filtered_records_df = filtered_records_df.applymap(lambda x: str(x).strip("[]'"))

# Save the filtered records to a new TSV file
filtered_records_df.to_csv(f'{args.dmb}_metadata.tsv', sep='\t', index=False)