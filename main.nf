nextflow.enable.dsl=2

include { gunzip_1 } from './modules/gunzip_1/module.nf'
include { gwas_vcf_regenie_1 } from './modules/gwas_vcf_regenie_1/module.nf'

workflow {
input1 = Channel.fromPath(params.gunzip_1.any_file)
input2 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_king_reference_data).splitCsv()
input3 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_pheno)
gunzip_1(input1)
gwas_vcf_regenie_1(gunzip_1.out.output1, input2, gunzip_1.out.output1, gunzip_1.out.output1, gunzip_1.out.output1, gunzip_1.out.output1, input3)
gwas_vcf_regenie_1(gunzip_1.out.output1, input2, gunzip_1.out.output1, gunzip_1.out.output1, gunzip_1.out.output1, gunzip_1.out.output1, input3)
}
