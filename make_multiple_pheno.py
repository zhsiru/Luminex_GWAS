import csv
import pandas as pd

df = pd.read_csv("Lum_case_inf30_max_norm.csv")
protein = list(df.columns.values)

outpath='/scratch/richards/sirui.zhou/GWAS_proteome/Luminex/'
for protein in df:
	x=pd.DataFrame(df, columns=["GWAS_ID","GWAS_ID",protein])
	outfile=outpath + protein + '_inf30_max_EUR' + '.pheno'
	x.to_csv(outfile, index=False, header=False, sep=' ')