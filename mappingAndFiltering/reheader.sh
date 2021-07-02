#! /bin/bash

module load bcftools


### SNP files
bcftools reheader -s reheader.Dmel.txt -o MS_Dmel.VarScan.newHeader.vcf MS_Dmel.VarScan.vcf
bcftools reheader -s reheader.Dsim.txt -o MS_Dsim.VarScan.newHeader.vcf MS_Dsim.VarScan.vcf
bcftools reheader -s reheader.Dsuz.txt -o MS_Dsuz.VarScan.newHeader.vcf MS_Dsuz.VarScan.vcf
bcftools reheader -s reheader.Dhyd.txt -o MS_Dhyd.VarScan.newHeader.vcf MS_Dhyd.VarScan.vcf
bcftools reheader -s reheader.Zind.txt -o MS_Zind.VarScan.newHeader.vcf MS_Zind.VarScan.vcf

### indel files 
bcftools reheader -s reheader.Dmel.txt -o MS_Dmel.VarScan.indel.newHeader.vcf MS_Dmel.VarScan.indel.vcf
bcftools reheader -s reheader.Dsim.txt -o MS_Dsim.VarScan.indel.newHeader.vcf MS_Dsim.VarScan.indel.vcf
bcftools reheader -s reheader.Dsuz.txt -o MS_Dsuz.VarScan.indel.newHeader.vcf MS_Dsuz.VarScan.indel.vcf
bcftools reheader -s reheader.Dhyd.txt -o MS_Dhyd.VarScan.indel.newHeader.vcf MS_Dhyd.VarScan.indel.vcf
bcftools reheader -s reheader.Zind.txt -o MS_Zind.VarScan.indel.newHeader.vcf MS_Zind.VarScan.indel.vcf
