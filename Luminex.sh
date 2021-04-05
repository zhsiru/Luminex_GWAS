
#!/bin/bash
#PBS -N Luminex
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=10G
#PBS -l vmem=10G
#PBS -l epilogue=/scratch/richards/sirui.zhou/GWAS_proteome/epilogue.sh

cd /scratch/richards/sirui.zhou/GWAS_proteome/
protein=($(awk '{ print $1}' Luminex_list.txt)) ###list of Luminex proteins from header of Lum_case_inf30_max_norm.csv
for ((i=0;i<${#protein[@]};++i)); do gcta64 --mbfile BQC_filelist.txt --fastGWA-lr --grm EUR/Lum_case_inf_max_norm_EUR --pheno Luminex/${protein[i]}_inf30_max_EUR.pheno --covar Lum_case_inf30_cov_EUR.txt --maf 0.05 --qcovar Lum_case_inf30_qcov_EUR.txt --thread-num 10 --out Luminex/${protein[i]}_inf30_max_EUR_lr; done

