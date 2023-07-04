## Raphael Eisenhofer 2023

import argparse
import csv
import pandas as pd
import requests
import json

#Load in CSV file
parser = argparse.ArgumentParser()
parser.add_argument('--dmb', required=True, help='DMB number')
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

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"{{DM_batch}} = '{args.dmb}'",

}

# Make the request to get the records
response = requests.get(url, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'dereped_mags.csv'

with open(output_file_path, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',')
    writer.writerow(['genome', 'completeness', 'contamination'])

    # Loop through the records and write each one to the csv file
    for record in records:
        # Get the ID of the record
        record_id = record['id']

        # Make a request to get the record data
        record_response = requests.get(f"{url}/{record_id}", headers=headers)

        # Extract the values of the completeness and contamination fields
        mag_name = record_response.json()['fields'].get('mag_name', '')
        completeness = record_response.json()['fields'].get('completeness', '')
        contamination = record_response.json()['fields'].get('contamination', '')

        # Write the row to the csv file
        row = [mag_name, completeness, contamination]
        writer.writerow(row)