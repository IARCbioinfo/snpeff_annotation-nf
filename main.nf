#!/usr/bin/env nextflow

params.input_folder_with_VCF_files = "${baseDir}/VCFs/"
params.reference_genome = "GRCh37.75"
params.dbNSF_path = "${baseDir}/dbNSFP4.1a.txt.gz"
params.dbSNP_path = "${baseDir}/dbsnp150.vcf.gz"
params.output_path = "${baseDir}/output"

process FilterInputFiles {
    tag "Sample ${sample}"

    input:
    path vcf_file

    output:
    val sample
    path("${sample}_filtered.vcf.gz")
    path("${sample}_filtered.vcf.gz.tbi")

    script:
    sample = vcf_file.toString().substring(0,vcf_file.toString().indexOf('.'))
    """
    # filter input VCF files by max number of alleles using plink
    plink2 --vcf $vcf_file \\
        --allow-extra-chr \\
        --max-alleles 2 \\
        --export vcf-4.2 \\
        --out ${sample}
    bcftools view -i "%FILTER='PASS'" ${sample}.vcf -Oz -o ${sample}_filtered.vcf.gz
    bcftools index --tbi ${sample}_filtered.vcf.gz
    """
}

process AnnotateWithRSID {
    tag "Sample $sample"

    input:
    val sample
    path vcf_file
    path vcf_file_tbi

    output:
    val sample
    path "${sample}_rsid.vcf.gz"
    path "${sample}_rsid.vcf.gz.tbi"

    script:
    def avail_mem = 4
    if (!task.memory) {
        log.info '[FullyAnnotateWithDbSNP] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    SnpSift -Xmx${avail_mem}g annotate -id -db ${params.dbSNP_path} ${vcf_file} > ${sample}_rsid.vcf
    bgzip ${sample}_rsid.vcf
    bcftools index --tbi ${sample}_rsid.vcf.gz
    """
}

process AnnotateWithImpact {
    tag "Sample $sample"
    
    input:
    val sample
    path vcf_file
    path vcf_file_tbi

    output:
    val sample
    path "${sample}_impact.vcf.gz"
    path "${sample}_impact.vcf.gz.tbi"

    script:
    def avail_mem = 4
    if (!task.memory) {
        log.info '[FullyAnnotateWithDbSNP] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    snpEff -Xmx${avail_mem}g ${params.reference_genome} ${vcf_file} > ${sample}_impact.vcf
    bgzip ${sample}_impact.vcf
    bcftools index --tbi ${sample}_impact.vcf.gz
    """
}

process FullyAnnotateWithDbSNP {
    tag "Sample $sample"
    
    input:
    val sample
    path vcf_file
    path vcf_file_tbi

    output:
    val sample
    path "${sample}_full.vcf.gz"
    path "${sample}_full.vcf.gz.tbi"

    script:
    def avail_mem = 4
    if (!task.memory) {
        log.info '[FullyAnnotateWithDbSNP] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    SnpSift \\
        -Xmx${avail_mem}g \\
        dbnsfp -f \\
        SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,gnomAD_genomes_AC,gnomAD_genomes_AF,gnomAD_genomes_AFR_AC,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMR_AC,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AC,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AC,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AC,gnomAD_genomes_FIN_AF,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_AMI_AC,gnomAD_genomes_AMI_AF,gnomAD_genomes_SAS_AC,gnomAD_genomes_SAS_AF,clinvar_id,clinvar_clnsig,clinvar_trait,REVEL_score,REVEL_rankscore \\
        -db ${params.dbNSF_path} ${vcf_file} > ${sample}_full.vcf
    bgzip ${sample}_full.vcf
    bcftools index --tbi ${sample}_full.vcf.gz
    """
}

process ExtractFields {
    tag "Sample $sample"

    input:
    val sample
    path vcf_file
    path vcf_file_tbi

    output:
    path "Full_annotation_${sample}.txt"

    script:
    def avail_mem = 4
    if (!task.memory) {
        log.info '[FullyAnnotateWithDbSNP] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    SnpSift \\
        -Xmx${avail_mem}g \\
        extractFields -s "," -e "." \\
        ${vcf_file} \\
        "CHROM" "POS" "ID" "REF" "ALT" \\
        "ANN[0].GENE" "ANN[0].IMPACT" "ANN[0].EFFECT" "GEN[0].GT" \\
        "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "GEN[*].GT" \\
        "dbNSFP_SIFT_pred[*]" "dbNSFP_Polyphen2_HDIV_pred[*]" \\
        "dbNSFP_Polyphen2_HVAR_pred[*]" "dbNSFP_LRT_pred[*]" \\
        "dbNSFP_gnomAD_genomes_AC[*]" "dbNSFP_gnomAD_genomes_AF[*]" \\
        "dbNSFP_gnomAD_genomes_AFR_AC[*]" "dbNSFP_gnomAD_genomes_AFR_AF[*]" \\
        "dbNSFP_gnomAD_genomes_AMR_AC[*]" "dbNSFP_gnomAD_genomes_AMR_AF[*]" \\
        "dbNSFP_gnomAD_genomes_ASJ_AC[*]" "dbNSFP_gnomAD_genomes_ASJ_AF[*]" \\
        "dbNSFP_gnomAD_genomes_EAS_AC[*]" "dbNSFP_gnomAD_genomes_EAS_AF[*]" \\
        "dbNSFP_gnomAD_genomes_FIN_AC[*]" "dbNSFP_gnomAD_genomes_FIN_AF[*]" \\
        "dbNSFP_gnomAD_genomes_NFE_AC[*]" "dbNSFP_gnomAD_genomes_NFE_AF[*]" \\
        "dbNSFP_gnomAD_genomes_AMI_AC[*]" "dbNSFP_gnomAD_genomes_AMI_AF[*]" \\
        "dbNSFP_gnomAD_genomes_SAS_AC[*]" "dbNSFP_gnomAD_genomes_SAS_AF[*]" \\
        "dbNSFP_REVEL_score[*]" "dbNSFP_REVEL_rankscore[*]" \\
        "dbNSFP_clinvar_id[*]" "dbNSFP_clinvar_clnsig[*]" \\
        "dbNSFP_clinvar_trait[*]" \\
        | awk '{print "'${sample}'\\t" \$0}' \\
        | sed -e '1s/${sample}/SAMPLE/' \\
        >> Full_annotation_${sample}.txt
    """
}

workflow {
    // Grab input VCF files
    file_channel = Channel.fromPath( params.input_folder_with_VCF_files + '/*vcf.gz', checkIfExists: true )
    // Launch the pipeline
    FilterInputFiles(file_channel) \
        | AnnotateWithRSID \
        | AnnotateWithImpact \
        | FullyAnnotateWithDbSNP \
        | ExtractFields \
        | collectFile(name: 'full_annotation.txt', \
            newLine: false, \
            keepHeader: true, \
            skip: 1, \
            sort: { file -> file.baseName }, \
            storeDir: params.output_path)
}
