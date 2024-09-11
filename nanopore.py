#!/usr/bin/env python

import os
import pandas as pd
import concurrent.futures

workflow_path = "/hps/nobackup/cochrane/ena/users/khadim/project/ena-covid-systematic-analysis-workflow"

params = {
    "SARS2_FA": os.path.join(workflow_path, "data/NC_045512.2.fa"),
    "SARS2_FA_FAI": os.path.join(workflow_path, "data/NC_045512.2.fa.fai"),
    "SECRETS": os.path.join(workflow_path, "prepro/projects_accounts.csv"),
    "INDEX": os.path.join(workflow_path, "prepro/nanopore.index.tsv"),
    "STOREDIR": os.path.join(workflow_path, "prepro/storeDir"),
    "OUTDIR": os.path.join(workflow_path, "prepro/results"),
    "STUDY": 'PRJEB45555'
}

def map_to_reference(run_accession, sample_accession, input_file, sars2_fasta, sars2_fasta_fai, projects_accounts_csv, study_accession, num_threads, STOREDIR, OUTDIR, workflow_path):
    # Read the line containing FTP credentials from the projects_accounts_csv file
    line = os.popen(f"grep {study_accession} {projects_accounts_csv}").read()
    ftp_id = line.split(',')[2]
    ftp_password = line.split(',')[4]
    
    # Download FASTQ file using wget
    if ftp_id == 'public':
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_1.fastq.gz {input_file}")
    else:
        os.system(f"cd {STOREDIR}; wget -t 0 -O {run_accession}_1.fastq.gz {input_file} --user={ftp_id} --password={ftp_password}")
    
    # Trim reads using cutadapt
    os.system(f"cd {STOREDIR}; cutadapt -u 30 -u -30 -o {run_accession}.trimmed.fastq {run_accession}_1.fastq.gz -m 75 -j {num_threads} --quiet")
    
    # Map reads to the reference genome using minimap2 and convert to BAM format using samtools
    os.system(f"cd {STOREDIR}; minimap2 -Y -t {num_threads} -x map-ont -a {sars2_fasta} {run_accession}.trimmed.fastq | samtools view -bF 4 - | samtools sort -@ {num_threads} - > {run_accession}.bam")
    os.system(f"cd {STOREDIR}; samtools index -@ {num_threads} {run_accession}.bam")
    
    # Generate pileup file
    os.system(f"cd {STOREDIR}; samtools mpileup -a -A -Q 30 -d 8000 -f {sars2_fasta} {run_accession}.bam > {run_accession}.pileup")
    
    # Generate coverage file
    os.system(f"cd {STOREDIR}; cat {run_accession}.pileup | awk '{{print $2,\",\"$3,\",\"$4}}' > {run_accession}.coverage")
    os.system(f"cd {STOREDIR}; bgzip -f  {run_accession}.coverage")

    # Index the BAM file
    os.system(f"cd {STOREDIR}; samtools index {run_accession}.bam")

    # Perform indel quality calibration using LoFreq
    os.system(f"cd {STOREDIR}; lofreq indelqual --dindel {run_accession}.bam -f {sars2_fasta} -o {run_accession}_fixed.bam")
    os.system(f"cd {STOREDIR}; samtools index {run_accession}_fixed.bam")
    
    # Convert BAM to VCF format using bam_to_vcf.py
    os.system(f"{workflow_path}/bin/bam_to_vcf.py -b {STOREDIR}/{run_accession}.bam -r {sars2_fasta} --mindepth 30 --minAF 0.1 -c {num_threads} -o {STOREDIR}/{run_accession}.vcf")
    
    # Filter VCF file using filtervcf.py and compress using bgzip -f 
    os.system(f"{workflow_path}/bin/filtervcf.py -i {STOREDIR}/{run_accession}.vcf -o {STOREDIR}/{run_accession}_filtered.vcf")
    os.system(f"cd {STOREDIR}; bgzip -f  {run_accession}.vcf")
    os.system(f"cd {STOREDIR}; bgzip -f  {run_accession}_filtered.vcf")
    
    # Index the filtered VCF file
    os.system(f"cd {STOREDIR}; tabix {run_accession}_filtered.vcf.gz")

    # Generate statistics for VCF files
    os.system(f"cd {STOREDIR}; bcftools stats {run_accession}_filtered.vcf.gz > {run_accession}.stat")
    
    # Convert VCF to consensus sequence using vcf2consensus.py and compress using bgzip -f 
    os.system(f"{workflow_path}/bin/vcf2consensus.py -v {STOREDIR}/{run_accession}_filtered.vcf.gz -d {STOREDIR}/{run_accession}.coverage -r {sars2_fasta} -o {STOREDIR}/{run_accession}_headless_consensus.fasta -dp 30 -n {run_accession}")
    os.system(f"{workflow_path}/bin/fix_consensus_header.py {STOREDIR}/{run_accession}_headless_consensus.fasta > {STOREDIR}/{run_accession}_consensus.fasta")
    os.system(f"cd {STOREDIR}; bgzip -f  {run_accession}_consensus.fasta")
    
    # Perform SNP annotation using SnpEff and compress the annotated VCF file
    os.system(f"cd {STOREDIR}; snpEff  -q -no-downstream -no-upstream -noStats NC_045512.2 {run_accession}.vcf > {run_accession}.annot.vcf")
    os.system(f"cd {STOREDIR}; bgzip -f  {run_accession}.annot.vcf")
    
    # Prepare output directory and move files
    os.system(f"cd {STOREDIR}; mkdir -p {run_accession}_output")
    os.system(f"cd {STOREDIR}; mv {run_accession}.annot.vcf {run_accession}.bam {run_accession}.coverage {run_accession}.stat {run_accession}.vcf.gz {run_accession}_output")
    
    # Create output tarball
    os.system(f"cd {STOREDIR}; tar -zcvf {run_accession}_output.tar.gz {run_accession}_output")

    # Create directory and copy files
    os.system(f"mkdir -p {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")
    os.system(f"cp {workflow_path}/bin/config.yaml {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")
    os.system(f"mv {STOREDIR}/{run_accession}_output.tar.gz {STOREDIR}/{run_accession}_filtered.vcf.gz {STOREDIR}/{run_accession}_consensus.fasta.gz {workflow_path}/prepro/results/{study_accession}/{run_accession}_output")

    # Cleanup
    os.system(f"cd {OUTDIR}; rm -rf {run_accession}*")
    


def ena_analysis_submit(run_accession, sample_accession, projects_accounts_csv, study_accession, workflowPath):
    # Read the line containing Webin ID and password from the projects_accounts_csv file
    line = os.popen(f"grep {study_accession} {projects_accounts_csv}").read()
    webin_id = line.split(',')[3]
    webin_password = line.split(',')[4]

    # Build paths for the output files
    output_tgz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_output.tar.gz"
    filtered_vcf_gz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_filtered.vcf.gz"
    consensus_fasta_gz_path = f"{workflowPath}/prepro/results/{study_accession}/{run_accession}_output/{run_accession}_consensus.fasta.gz"

    # Submit files using ena-analysis-submitter.py
    if study_accession == 'PRJEB45555':
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB43947 -s {sample_accession} -r {run_accession} -f {output_tgz_path} -a PATHOGEN_ANALYSIS -au {webin_id} -ap {webin_password} > {output_tgz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB45554 -s {sample_accession} -r {run_accession} -f {filtered_vcf_gz_path} -a COVID19_FILTERED_VCF -au {webin_id} -ap {webin_password} > {filtered_vcf_gz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p PRJEB45619 -s {sample_accession} -r {run_accession} -f {consensus_fasta_gz_path} -a COVID19_CONSENSUS -au {webin_id} -ap {webin_password} > {consensus_fasta_gz_path}_submit.txt")
    else:
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {output_tgz_path} -a PATHOGEN_ANALYSIS -au {webin_id} -ap {webin_password} > {output_tgz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {filtered_vcf_gz_path} -a COVID19_FILTERED_VCF -au {webin_id} -ap {webin_password} > {filtered_vcf_gz_path}_submit.txt")
        os.system(f"{workflowPath}/bin/ena-analysis-submitter.py -p {study_accession} -s {sample_accession} -r {run_accession} -f {consensus_fasta_gz_path} -a COVID19_CONSENSUS -au {webin_id} -ap {webin_password} > {consensus_fasta_gz_path}_submit.txt")

def process_data(row):
    try:
        if ';' in row['fastq_ftp']:
            fastq_urls = row['fastq_ftp'].split(';')[0]
        else:
            fastq_urls = row['fastq_ftp']
        study_accession = params['STUDY']
        run_accession = row['run_accession']
        sample_accession = row['sample_accession']
        output_folder = f"{workflow_path}/prepro/results/{study_accession}/{run_accession}_output"
        
        with open('accession_submitted.txt', 'r') as file:
            accession_submitted = file.read()

        # If output files don't exist, run map_to_reference
        if run_accession not in accession_submitted :
            map_to_reference(row['run_accession'], row['sample_accession'], fastq_urls, params['SARS2_FA'], params['SARS2_FA_FAI'], params['SECRETS'], params['STUDY'], num_threads=8, STOREDIR=params['STOREDIR'], OUTDIR=params['OUTDIR'], workflow_path=workflow_path)
            ena_analysis_submit(row['run_accession'], row['sample_accession'], params['SECRETS'], params['STUDY'], workflow_path)
            os.system(f'cat {output_folder}/*txt > {output_folder}/submit.txt')
            with open(f'{output_folder}/submit.txt', 'r') as submit_file:
                submit = submit_file.read()
            if 'error' not in submit.lower():
                with open('accession_submitted.txt', 'a') as f:
                    f.write(f'{run_accession}\n')

            os.system(f'rm -rf {output_folder}')
    except Exception as e:
        os.system(f'echo {e} >> getion_error_ont.err')

def main():
    index = pd.read_csv(params['INDEX'], sep='\t')
    
    # Number of tasks to run in parallel
    max_workers = 50
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_data, row) for _, row in index.iterrows()]
        
        # Wait for all tasks to finish
        for future in concurrent.futures.as_completed(futures):
            pass  # Can be extended for handling results or exceptions

if __name__ == "__main__":
    main()
