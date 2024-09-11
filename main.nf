#!/usr/bin/env nextflow

//params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
//params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"
//params.SECRETS = "gs://prj-int-dev-covid19-nf-gls/prepro/projects_accounts.csv"

params.INDEX = "/hps/nobackup/cochrane/ena/users/azyoud/software/systematicAnalysisWorkflow/covid-sequence-analysis-workflow/prepro/illumina.index.tsv"
params.STOREDIR = "/hps/nobackup/cochrane/ena/users/azyoud/software/systematicAnalysisWorkflow/covid-sequence-analysis-workflow/prepro/storeDir"
params.OUTDIR = "/hps/nobackup/cochrane/ena/users/azyoud/software/systematicAnalysisWorkflow/covid-sequence-analysis-workflow/prepro/results"

params.STUDY = 'PRJEB45555'

nextflow.enable.dsl = 2

process map_to_reference {
    storeDir params.STOREDIR

    cpus 6
    memory '8 GB'
    container 'davidyuyuan/ena-sars-cov2-illumina:1.0'

    input:
    tuple val(run_accession), val(sample_accession), file(input_file_1), file(input_file_2) //from samples_ch
    path(sars2_fasta)
    path(sars2_fasta_fai)
    path(projects_accounts_csv)
    val(study_accession)

    output:
    val(run_accession)
    val(sample_accession)
    file("${run_accession}_output.tar.gz")
    file("${run_accession}_filtered.vcf.gz")
    file("${run_accession}_consensus.fasta.gz")

    script:
    """
    line="\$(grep ${study_accession} ${projects_accounts_csv})"
    ftp_id="\$(echo \${line} | cut -d ',' -f 3)"
    ftp_password="\$(echo \${line} | cut -d ',' -f 6)"
    
    if [ "\${ftp_id}" = 'public' ]; then
        wget -t 0 -O ${run_accession}_1.fastq.gz \$(cat ${input_file_1})
        wget -t 0 -O ${run_accession}_2.fastq.gz \$(cat ${input_file_2})
    else
        wget -t 0 -O ${run_accession}_1.fastq.gz \$(cat ${input_file_1}) --user=\${ftp_id} --password=\${ftp_password}
        wget -t 0 -O ${run_accession}_2.fastq.gz \$(cat ${input_file_2}) --user=\${ftp_id} --password=\${ftp_password}
    fi
    trimmomatic PE ${run_accession}_1.fastq.gz ${run_accession}_2.fastq.gz ${run_accession}_trim_1.fq \
    ${run_accession}_trim_1_un.fq ${run_accession}_trim_2.fq ${run_accession}_trim_2_un.fq \
    -summary ${run_accession}_trim_summary -threads ${task.cpus} \
    SLIDINGWINDOW:5:30 MINLEN:50

    bwa index ${sars2_fasta}
    bwa mem -t ${task.cpus} ${sars2_fasta} ${run_accession}_trim_1.fq ${run_accession}_trim_2.fq | samtools view -bF 4 - | samtools sort - > ${run_accession}_paired.bam
    bwa mem -t ${task.cpus} ${sars2_fasta} <(cat ${run_accession}_trim_1_un.fq ${run_accession}_trim_2_un.fq) | samtools view -bF 4 - | samtools sort - > ${run_accession}_unpaired.bam
    samtools merge ${run_accession}.bam ${run_accession}_paired.bam ${run_accession}_unpaired.bam
    rm ${run_accession}_paired.bam ${run_accession}_unpaired.bam

    samtools mpileup -a -A -Q 30 -d 8000 -f ${sars2_fasta} ${run_accession}.bam > ${run_accession}.pileup
    cat ${run_accession}.pileup | awk '{print \$2,","\$3,","\$4}' > ${run_accession}.coverage

    samtools index ${run_accession}.bam
    lofreq indelqual --dindel ${run_accession}.bam -f ${sars2_fasta} -o ${run_accession}_fixed.bam
    samtools index ${run_accession}_fixed.bam
    lofreq call-parallel --no-default-filter --call-indels --pp-threads ${task.cpus} -f ${sars2_fasta} -o ${run_accession}.vcf ${run_accession}_fixed.bam
    lofreq filter --af-min 0.25 -i ${run_accession}.vcf -o ${run_accession}_filtered.vcf
    bgzip ${run_accession}.vcf
    bgzip ${run_accession}_filtered.vcf
    tabix ${run_accession}.vcf.gz
    bcftools stats ${run_accession}.vcf.gz > ${run_accession}.stat

    snpEff -q -no-downstream -no-upstream -noStats NC_045512.2 ${run_accession}.vcf > ${run_accession}.annot.vcf
    # vcf_to_consensus.py -dp 10 -af 0.25 -v ${run_accession}.vcf.gz -d ${run_accession}.coverage -o ${run_accession}_consensus.fasta -n ${run_accession} -r ${sars2_fasta}
    vcf_to_consensus.py -dp 10 -af 0.25 -v ${run_accession}.vcf.gz -d ${run_accession}.coverage -o headless_consensus.fasta -n ${run_accession} -r ${sars2_fasta}
    fix_consensus_header.py headless_consensus.fasta > ${run_accession}_consensus.fasta
    bgzip ${run_accession}_consensus.fasta

    mkdir -p ${run_accession}_output
    mv ${run_accession}_trim_summary ${run_accession}.annot.vcf ${run_accession}.bam ${run_accession}.coverage ${run_accession}.stat ${run_accession}.vcf.gz ${run_accession}_output
    tar -zcvf ${run_accession}_output.tar.gz ${run_accession}_output
    """
}

// Same process duplicated in both pipelines due to a limitation in Nextflow
process ena_analysis_submit {
    publishDir params.OUTDIR, mode: 'move'
    storeDir params.STOREDIR

    container 'davidyuyuan/ena-analysis-submitter:2.0'

    input:
    val(run_accession)
    val(sample_accession)
    file(output_tgz)
    file(filtered_vcf_gz)
    file(consensus_fasta_gz)
    path(projects_accounts_csv)
    val(study_accession)

    output:
    file("${run_accession}_output/${study_accession}/${run_accession}_output.tar.gz")
    file("${run_accession}_output/${study_accession}/${run_accession}_filtered.vcf.gz")
    file("${run_accession}_output/${study_accession}/${run_accession}_consensus.fasta.gz")
    file("${run_accession}_output/${study_accession}/successful_submissions.txt")

    script:
    """
    line="\$(grep ${study_accession} ${projects_accounts_csv})"
    webin_id="\$(echo \${line} | cut -d ',' -f 4)"
    webin_password="\$(echo \${line} | cut -d ',' -f 5)"
    
    mkdir -p ${run_accession}_output/${study_accession}
    cp /hps/nobackup/cochrane/ena/users/azyoud/software/systematicAnalysisWorkflow/covid-sequence-analysis-workflow/bin/config.yaml ${run_accession}_output/${study_accession}
    if [ "${study_accession}" = 'PRJEB45555' ]; then
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p PRJEB43947 -s ${sample_accession} -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p PRJEB45554 -s ${sample_accession} -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p PRJEB45619 -s ${sample_accession} -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password}
    else
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t true -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password}
    fi
    mv ${output_tgz} ${filtered_vcf_gz} ${consensus_fasta_gz} ${run_accession}_output/${study_accession}
    """
}

workflow {
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp.split(';')[0], 'ftp://' + row.fastq_ftp.split(';')[1]) }

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI, params.SECRETS, params.STUDY)
    ena_analysis_submit(map_to_reference.out, params.SECRETS, params.STUDY)
}
