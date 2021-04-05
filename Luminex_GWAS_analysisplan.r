###1. Prepare phenotype/covariate data for fastGWA###

library(dplyr)  
library(data.table)
library(ggplot2)
library(reshape2)
library(SomaDataIO)

###"Luminex_0401.txt" is a masterfile with clinical information, Luminex data and genotype ID###
Lum <- read.table("Luminex_0401.txt",header=T,fill=TRUE)
Lum <- Lum[order(Lum$Date_Sampling), ]
###take only COVID-19 case with Luminex data
Lum_case <- Lum[ which(Lum$GWAS_ID != "NA" 
                       & Lum$CCL2 != "NA" 
                       & Lum$COVID.19_Status == "POSITIF"), ]
###take only samples during infection (Day of Symptom onset (DSO)<30)					   
Lum$DSO <- as.numeric(Lum$DSO)
Lum_case_inf <- Lum_case[ which(Lum_case$DSO < 30), ]
###select max value sample per individual and log+standardize transform
Lum_case_inf_max <- Lum_case_inf %>%
  group_by(GWAS_ID) %>%  
  summarise_at(vars(30:57), max)
Lum_case_inf_max_norm <- Lum_case_inf_max
Lum_case_inf_max_norm[c(2:29)] <- scale(log(Lum_case_inf_max[c(2:29)]))

###If want to use baseline sample###
Lum_case_inf_1st <- Lum_case_inf %>% 
  group_by(GWAS_ID) %>% 
  filter(VAP == first(VAP))

###take a look of distribution after standardize
ggplot(melt(Lum_case_inf_max_norm),aes(x=value)) + geom_histogram() + facet_wrap(~variable)

###export phenotype and sample file
write.csv(Lum_case_inf_max_norm, file="Lum_case_inf30_max_norm.csv",row.names=F)
Lum_case_inf_max_norm_ind <- Lum_case_inf_max_norm[c(1,1)]
write.table(Lum_case_inf_max_norm_ind, file="Lum_case_inf_max_norm.ind",col.names=F,row.names=F,quote=F,sep=" ")

###2. run gcta to create grm and PC, we already have a European only grm file created. 
###gcta64  --bfile EUR/grm_final_EUR  --keep Lum_case_inf_max_norm.ind  --autosome  --maf 0.01 --make-grm  --out EUR/Lum_case_inf_max_norm_EUR
###gcta64  --grm EUR/Lum_case_inf_max_norm_EUR --pca 5 --out EUR/Lum_case_inf_max_norm_EUR_PC

###3. Create covariate and qcovariate file for fastGWA
EUR_Lum_PC <- fread("Lum_case_inf_max_norm_EUR_PC.eigenvec", header=F, select = c(2:7))
colnames(EUR_Lum_PC) <- c("GWAS_ID","PC1","PC2","PC3","PC4","PC5")
###take ID, sex and center and remove duplicates
Lum_case_inf30_cov <- unique(Lum_case_inf[c(3,4,8)])
###merge and format
Lum_case_inf30_cov_EUR <- merge(Lum_case_inf30_cov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1,2,3)]
###take ID, age and remove duplicates
Lum_case_inf30_qcov <- unique(Lum_case_inf[c(3,7)])
###merge with PC and format
Lum_case_inf30_qcov_EUR <- merge(Lum_case_inf30_qcov, EUR_Lum_PC, by = "GWAS_ID")[c(1,1:7)]
###export covariate and qcovariate files
write.table(Lum_case_inf30_cov_EUR, file="Lum_case_inf30_cov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(Lum_case_inf30_qcov_EUR, file="Lum_case_inf30_qcov_EUR.txt",col.names=T,row.names=F,quote=F,sep="\t")


###4. run make_multiple_pheno.py to split phenotype file for fastGWA

###5. run fastGWA using Luminex.sh 
