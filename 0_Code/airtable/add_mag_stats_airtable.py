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
    # Set the cell data you want to update
    data = {
        'fields': {
            'mag_name': str(row['mag_name']),
            'eha_number': str(row['eha_number']),
            'GTDB_version': chr(row['GTDB_version']),
            'domain': str(row['domain']),
            'phylum': str(row['phylum']),
            'class': str(row['class']),
            'order': str(row['order']),
            'family': str(row['family']),
            'genus': str(row['genus']),
            'species': str(row['species']),
            'closest_placement_ani': int(row['closest_placement_ani']),
            'completeness': int(row['completeness']),
            'contamination': int(row['contamination']),
            'size': int(row['size']),
            'GC': int(row['GC']),
            'N50': int(row['N50']),
            'contigs': int(row['contigs'])
        }
    }

    # Send a POST request to create a new record
    response = requests.post(url, headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)
