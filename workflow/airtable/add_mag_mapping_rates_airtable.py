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
url = 'https://api.airtable.com/v0/appWbHBNLE6iAsMRV/tblWDyQmM9rQ9wq57'
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
            'PR_batch_static': str(row['PR_batch_static']),
            'EHI_sample_static': str(row['EHI_sample_static']),
            'DM_batch_static': str(row['DM_batch_static']),
            'MAG_mapping_percentage': float(row['MAG_mapping_percentage']),
        }
    }

    # Send a POST request to create a new record
    response = requests.post(url, headers=headers, data=json.dumps(data))

    # Print the response status code
    print(response.status_code)
