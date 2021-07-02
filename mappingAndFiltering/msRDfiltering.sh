#! /bin/bash

module load gcc/7.1.0
module load vcftools/0.1.15


vcftools --vcf MS_Dmel.VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf --positions MS_Dmel_RDkeep.txt --recode --recode-INFO-all --out MS_Dmel.VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt

vcftools --vcf MS_Dsim.VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf --positions MS_Dsim_RDkeep.txt --recode --recode-INFO-all --out MS_Dsim.VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt

vcftools --vcf MS_Dhyd.VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf --positions MS_Dhyd_RDkeep.txt --recode --recode-INFO-all --out MS_Dhyd.VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt

vcftools --vcf MS_Dsuz.VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf --positions MS_Dsuz_RDkeep.txt --recode --recode-INFO-all --out MS_Dsuz.VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt

vcftools --vcf MS_Zind.VarScan.newHeader.noIndelReg.noRep.biAllelic.recode.vcf --positions MS_Zind_RDkeep.txt --recode --recode-INFO-all --out MS_Zind.VarScan.newHeader.noIndelReg.noRep.biAllelic.RDfilt
