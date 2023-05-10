import argparse
import requests
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('--seb', required=True, help='SEB batch number')
args = parser.parse_args()

# Read the API key from config file
with open('/projects/ehi/data/.airtable_api_key.json') as f:
    config = json.load(f)
api_key = config['api_key']

# Set up the Airtable API endpoint
AIRTABLE_API_ENDPOINT = f'https://api.airtable.com/v0/appQpr6MxnaiVHsHy/tbldKeYbUhEYC9NGQ'

# Set up the request headers with the API key
headers = {
    'Authorization': f'Bearer {api_key}'
}

# Set up the query parameters to filter the records
query_params = {
    'filterByFormula': f"AND({{SE_batch}} = '{args.seb}', {{ena_specimen_exists}} = 'FALSE')",
    'view': 'viwQJJtpZknO0TM8L'
}

# Make the request to get the records
response = requests.get(AIRTABLE_API_ENDPOINT, params=query_params, headers=headers)

# Convert the response to a JSON object
data = response.json()

# Extract the records from the JSON object
records = data['records']

# Set up the output TSV file
output_file_path = 'output.tsv'

# Open the output file in write mode and set up the CSV writer
with open(output_file_path, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t')

    # Write the header row
    header_row = list(records[0]['fields'].keys())
    writer.writerow(header_row)

    # Write the data rows
    for record in records:
        # Extract the field values from the record
        fields = record['fields']

        # Create a list of values for the row
        row = list(fields.values())

        # Write the row to the output file
        writer.writerow(row)