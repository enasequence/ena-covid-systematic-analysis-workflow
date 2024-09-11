import argparse
import requests
import io
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta  # Import the datetime module

BASE_PORTAL_API_SEARCH_URL = 'https://www.ebi.ac.uk/ena/portal/api/search'
#FIELDS = ['run_accession', 'sample_accession', 'instrument_platform', 'instrument_model', 'fastq_aspera',
#          'fastq_bytes', 'fastq_ftp', 'fastq_galaxy', 'fastq_md5', 'first_created', 'first_public', 'country',
#          'collection_date', 'isolate', 'strain']

FIELDS = ['run_accession', 'sample_accession','fastq_ftp', 'first_created']



def get_args():
    parser = argparse.ArgumentParser(description='\033[93mThis script is used to download new ONT and Illumina submissions (runs) of COVID-19. \033[0m')
    parser.add_argument('-p', '--platform', type=str, help='Instrument platform', required=True)
    return parser.parse_args()


def get_url(instrument_platform):
    return f"{BASE_PORTAL_API_SEARCH_URL}?result=read_run&query=tax_tree(2697049)%20AND%20instrument_platform%3D%22{instrument_platform}%22&fields={'%2C'.join(FIELDS)}&format=tsv&limit=0"


def process_data(response_content):
    return pd.read_csv(io.StringIO(response_content.decode('UTF-8')), sep="\t", low_memory=False)

def clean_data(data):
    clean_data = data.dropna(subset=['run_accession', 'sample_accession', 'fastq_ftp'])
    return clean_data



def save_data(data, instrument_platform, date):
    if instrument_platform == 'oxford_nanopore':
        instrument_platform = 'nanopore'
    output_file = f"prepro/{instrument_platform}.index.{date}.tsv"
    cleaned_data = clean_data(data)
    cleaned_data.to_csv(output_file, sep="\t", index=False)
    os.system(f"cd prepro ; cat {instrument_platform}.index.{date}.tsv > {instrument_platform}.index.tsv")


def fetch_and_save_data(platform):
    # Calculate today's date minus two months
    date = (datetime.now() - timedelta(days=2)).strftime('%Y-%m-%d')

    url = get_url(platform)
    response = requests.get(url)
    if response.status_code == 200:
        data = process_data(response.content)
        data['first_created'] = pd.to_datetime(data['first_created'])
        filtered_data = data[data["first_created"] >= date]
        save_data(filtered_data, platform, date)
    else:
        print(f"Failed to fetch data for {platform} platform. Status code: {response.status_code}")


def main():
    args = get_args()
    fetch_and_save_data(args.platform)


if __name__ == '__main__':
    main()
