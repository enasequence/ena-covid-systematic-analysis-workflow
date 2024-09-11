#!/usr/bin/env python

__author__ = "Khadim GUEYE"

import argparse
import subprocess
import hashlib
import yaml , os , sys
from datetime import datetime
from lxml import etree


def main():
    args = get_args()
    
    # Read config.yaml file
    output_dir = os.path.dirname(args.file)
    config_path = os.path.join(output_dir, 'config.yaml')
    with open(config_path, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    # Generate XML using arguments and configuration
    sub_path = generate_submission_xml(args, config)
    xml_path = generate_analysis_xml(args, config)
    
    # file upload
    # upload_file(args)


    # xml file upload 
    submit_to_ENA(sub_path, xml_path, args)

def generate_submission_xml(args, config):
    output_dir = os.path.dirname(args.file)
    submission_filename = os.path.join(output_dir, "submission.xml")

    # Define the XML structure for submission
    #submission_root = etree.Element("SUBMISSION_SET")
    submission = etree.Element("SUBMISSION")
    actions = etree.SubElement(submission, "ACTIONS")
    action = etree.SubElement(actions, "ACTION")
    etree.SubElement(action, "ADD")

    # Create the XML tree
    submission_tree = etree.ElementTree(submission)

    # Write XML content to file
    with open(submission_filename, 'wb') as xml_sub_file:
        submission_tree.write(xml_sub_file, pretty_print=True, xml_declaration=True, encoding='UTF-8')

    return submission_filename

def generate_analysis_xml(args, config):
    output_dir = os.path.dirname(args.file)
    analysis_filename = os.path.join(output_dir, "analysis.xml")

    # Define the XML structure for analysis
    analysis_root = etree.Element("ANALYSIS_SET")
    analysis = etree.SubElement(analysis_root, "ANALYSIS")
    analysis.set("alias", f"{config['ALIAS']}_{args.project}_{args.sample}_{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}")
    analysis.set("center_name", config['CENTER_NAME'])
    analysis.set("analysis_date", datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))

    # Add elements to analysis XML
    etree.SubElement(analysis, "TITLE").text = config['TITLE']
    etree.SubElement(analysis, "DESCRIPTION").text = config['DESCRIPTION']
    etree.SubElement(analysis, "STUDY_REF").set("accession", args.project)
    etree.SubElement(analysis, "SAMPLE_REF").set("accession", args.sample)
    etree.SubElement(analysis, "RUN_REF").set("accession", args.run)
    analysis_type = etree.SubElement(analysis, "ANALYSIS_TYPE")
    etree.SubElement(analysis_type, args.analysis_type)

    files = etree.SubElement(analysis, "FILES")
    analysis_file = etree.SubElement(files, "FILE")
    filename = os.path.basename(args.file)
    analysis_file.set("filename",filename)
    if filename.endswith('.fasta.gz'):
        analysis_file.set("filetype", "fasta")
    elif filename.endswith('.vcf.gz'):
        analysis_file.set("filetype", "vcf")
    else:
        analysis_file.set("filetype", "other")
    analysis_file.set("checksum_method", "MD5")
    analysis_file.set("checksum", calculate_md5(args.file))  
    analysis_attributes = etree.SubElement(analysis, "ANALYSIS_ATTRIBUTES")
    for key, value in config.items():
        attribute = etree.SubElement(analysis_attributes, "ANALYSIS_ATTRIBUTE")
        etree.SubElement(attribute, "TAG").text = key
        etree.SubElement(attribute, "VALUE").text = value

    # Create the XML tree
    analysis_tree = etree.ElementTree(analysis_root)

    # Write XML content to file
    with open(analysis_filename, 'wb') as xml_file:
        analysis_tree.write(xml_file, pretty_print=True, xml_declaration=True, encoding='UTF-8')

    return analysis_filename


def calculate_md5(file_path):
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()

def upload_file(args):
	filename = os.path.basename(args.file)
	output_dir = os.path.dirname(args.file)
	curl_command= f'cd {output_dir} ; lftp -u {args.analysis_username},{args.analysis_password} -e "mput {filename}; bye" webin2.ebi.ac.uk'
	try:
		os.system(curl_command)
		#print(f"\033[93m{args.file} uploaded.\033[0m")
	except subprocess.CalledProcessError as e:
		print(f"\033[91mError: Submission failed.\033[0m\nException: {e}")
		sys.exit()

def submit_to_ENA(sub_path, xml_path, args):
	filename = os.path.basename(args.file)
	output_dir = os.path.dirname(args.file)
	curl_command_upload= f'cd {output_dir} ; lftp -u {args.analysis_username},{args.analysis_password} -e "mput {filename}; bye" webin2.ebi.ac.uk'
	# Construct the curl command
	if args.test:
        
		curl_command = [
		f'curl',
		'-u', f'{args.analysis_username}:{args.analysis_password}',
		'-F', f'SUBMISSION=@{sub_path}',
		'-F', f'ANALYSIS=@{xml_path}',
		'https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/'
		]
	else:
		curl_command = [
		f'curl',
		'-u', f'{args.analysis_username}:{args.analysis_password}',
		'-F', f'SUBMISSION=@{sub_path}',
		'-F', f'ANALYSIS=@{xml_path}',
		'https://www.ebi.ac.uk/ena/submit/drop-box/submit/'
		]

	# Execute the curl command
	try:
		os.system(curl_command_upload)
		subprocess.run(curl_command, check=True)
		#print("\033[93mSubmission successful.\033[0m")
	except subprocess.CalledProcessError as e:
		print(f"\033[91mError: Submission failed.\033[0m\nException: {e}")
		sys.exit()

def get_args():
    '''
    Define and obtain script arguments
    :return: Arguments object
    '''
    parser = argparse.ArgumentParser(prog='ena-analysis-submitter.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""
        + ============================================================ +
        |  European Nucleotide Archive (ENA) Analysis Submission Tool  |
        |                                                              |
        |  Tool to submit Pathogen analysis objects to an ENA project  |
        |  , mainly in the data hub context.                           |
        + =========================================================== +
        """)
    parser.add_argument('-p', '--project', help='Valid ENA project accession to submit analysis to (e.g. PRJXXXXXXX)', type=str, required=True)
    parser.add_argument('-s', '--sample', help='ENA sample accessions/s to link with the analysis submission', required=True)
    parser.add_argument('-r', '--run', help='ENA run accession/s to link with the analysis submission', required=True)
    parser.add_argument('-f', '--file', help='Files of analysis to submit to the project, accepts a list of files (e.g. path/to/file1.csv.gz,path/to/file2.txt.gz)', type=str, required=True)
    parser.add_argument('-a', '--analysis_type', help='Type of analysis to submit. Options: PATHOGEN_ANALYSIS, COVID19_CONSENSUS, COVID19_FILTERED_VCF, PHYLOGENY_ANALYSIS, FILTERED_VARIATION, SEQUENCE_CONSENSUS', choices=['PATHOGEN_ANALYSIS', 'COVID19_CONSENSUS', 'COVID19_FILTERED_VCF', 'PHYLOGENY_ANALYSIS', 'FILTERED_VARIATION', 'SEQUENCE_CONSENSUS'], required=True)         # Can add more options if you wish to share more analysis types
    parser.add_argument('-au', '--analysis_username', help='Valid Webin submission account ID (e.g. Webin-XXXXX) used to carry out the submission', type=str, required=True)
    parser.add_argument('-ap', '--analysis_password', help='Password for Webin submission account', type=str, required=True)
    parser.add_argument('-t', '--test', action='store_true', help='Specify whether to use ENA test server for submission.')
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()
