import argparse
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--report', required=True, help='Path to genome report file')
args = parser.parse_args()

#Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

#Set variables
url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblMzd3oyaJhdeQcs'
headers = {
    'Authorization': f'Bearer {api_key}',
    'Content-Type': 'application/json'
}

#Read in the TSV file using pandas
df = pd.read_csv(args.report, sep='\t')

# Loop through each row in the dataframe
for i, row in df.iterrows():
    params = {
        'maxRecords': 1
    }
    response = requests.get(url, headers=headers, params=params)
    data = response.json()
    record_id = data['records'][0]['id']

    # Set the cell data you want to update
    data = {
        'fields': {
            'mag_name': int(row['mag_name']),
            'eha_number': int(row['eha_number']),
            'GTDB_version': int(row['GTDB_version']),
            'domain': int(row['domain']),
            'phylum': int(row['phylum']),
            'class': int(row['class']),
            'order': int(row['order']),
            'family': int(row['family']),
            'genus': int(row['genus']),
            'species': int(row['species']),
            'closest_placement_ani': int(row['closest_placement_ani']),
            'completeness': int(row['completeness']),
            'contamination': int(row['contamination']),
            'size': int(row['size']),
            'GC': int(row['GC']),
            'N50': int(row['N50']),
            'contigs': int(row['contigs'])
        }
    }

    # Send a PATCH request to update the record
    response = requests.patch(f'{url}/{record_id}', headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)
