setwd("/Users/1Air/Documents/HS_630/BRCA_data")
#source("https://bioconductor.org/biocLite.R")
library(RTCGAToolbox)
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
library("TCGAbiolinks")
library(dplyr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("RTCGA")
library(RTCGA)
#install.packages("magrittr")
library(magrittr)

# import 3 datasets
# 1. Clinical data
#######################

downloadTCGA(cancerTypes = "BRCA",
             dataSet = "Merge_Clinical.Level_1",
             destDir = "BRCA_data",
             date = "2015-11-01")

readTCGA(path = file.path("BRCA_data",
                          grep("clinical_clin_format.txt",
                               list.files("BRCA_data/",
                                          recursive = TRUE),
                               value = TRUE)
),
dataType = "clinical") -> BRCA.clinical.20151101

###########################
#2. Microarray and RNAseq data

brca_data <- getFirehoseData(dataset = "BRCA", runDate = "20150821", 
                             gistic2_Date="20150821", RNAseq_Gene="TRUE", CNV_SNP = "TRUE", 
                             CNA_SNP = "TRUE", mRNA_Array ="TRUE")

###########################

# 3. my_data htseq rnaseq

dir.create('data2')

# downloading rnaseq data
downloadTCGA( cancerTypes = 'BRCA', 
              dataSet = 'rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level',
              destDir = 'data2' )

# shortening paths and directories
list.files( 'data2/') %>% 
  file.path( 'data2', .) %>%
  file.rename( to = substr(.,start=1,stop=50))

# reading data
list.files( 'data2/') %>% 
  file.path( 'data2', .) -> folder

folder %>%
  list.files %>%
  file.path( folder, .) %>%
  grep( pattern = 'illuminahiseq', x = ., value = TRUE) -> pathRNA
readTCGA( path = pathRNA, dataType = 'rnaseq' ) -> my_data

#########################
# start preprocessing
########################
# extract microarray and rnaseq data into large matrices

#brca_data.cnvsnp <- getData(brca_data, type = "CNVSNP")
#brca_data.cnasnp <- getData(brca_data, type = "CNASNP")
brca_data.mrnaArray <- getData(brca_data,type="mRNAArray",platform = 1)
brca_data.rna_seq <- getData(brca_data,type="RNASeqGene")

# extract clinical data
#clinical_data2 <- data.frame(brca_data@Clinical$years_to_birth, brca_data@Clinical$pathologic_stage, brca_data@Clinical$gender, brca_data@Clinical$histological_type, brca_data@Clinical$race, brca_data@Clinical$ethnicity, brca_data@RNASeqGene)
clinical_data <- data.frame(BRCA.clinical.20151101$patient.bcr_patient_barcode, BRCA.clinical.20151101$patient.breast_carcinoma_progesterone_receptor_status, BRCA.clinical.20151101$patient.breast_carcinoma_estrogen_receptor_status, BRCA.clinical.20151101$patient.menopause_status, BRCA.clinical.20151101$patient.histological_type, BRCA.clinical.20151101$patient.race, BRCA.clinical.20151101$patient.age_at_initial_pathologic_diagnosis, stringsAsFactors = FALSE)
# preprocess clinical data
colnames(clinical_data) <- c("barcode", "progestrone_status", "estrogen_status", "menopause_status", "histo_type", "race", "age")

clinical_data$age <- as.numeric(clinical_data$age)
clinical_data$barcode <- as.character(clinical_data$barcode)
clinical_data$barcode_short <- substr(clinical_data$barcode, 1, 12)
cytokines_long_bcr <- as.character(brca_cytokines$bcr_patient_barcode)
brca_cytokines$barcode_short <- substr(brca_cytokines$bcr_patient_barcode, 1, 12)

#########
# subset Estrogen receptor positive and negative to extract barcodes

er_pos <- subset(clinical_data, clinical_data$estrogen_status == 'positive' )
er_neg <- subset(clinical_data, clinical_data$estrogen_status == 'negative')
er_pos_bcr <- c(er_pos$barcode)
er_neg_bcr <- c(er_neg$barcode)
er_pos_bcr <- toupper(er_pos_bcr)
er_neg_bcr <- toupper(er_neg_bcr)

#####################################
# extract select cytokines gene expression quantification

# larger my_data rna dataset on 1212 obs
brca_cytokines <- my_data[, c("bcr_patient_barcode", "IL1B|3553","IL2|3558","IL2RA|3559","IL4|3565","IL6|3569","IL8|3576", 
                              "IL10|3586","IL12A|3592","IL12B|3593","IL12RB1|3594","IL12RB2|3595","IL13|3596","IL13RA1|3597","IL13RA2|3598",
                              "IL18|3606","IFNA10|3446","IFNA13|3447","IFNA14|3448","IFNA16|3449","IFNA17|3451","IFNA1|3439",
                              "IFNB1|3456","IFNG|3458","IRF1|3659","IRF2BP1|26145","IRF2BP2|359948","IRF2|3660",
                              "TGFBR3|7049","TGFBR2|7048","TGFBI|7045","TGFB3|7043","TGFB2|7042","TGFB1|7040"
)]


# microarray cytokines gene expression on 1098 obs

cytokines_mrna <- brca_data.mrnaArray[c("IL1B","IL4","IL2","IL2RA","IL4","IL6","IL8","IL10",
"IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL13RA1","IL13RA2","IL18","IFNA1","IFNB1",
"IFNG","IRF1","IRF2"
),]

# rna cytokines gene expression on 1098 obs

cytokines_rna <- brca_data.rna_seq[c("IL1B","IL4","IL2","IL2RA","IL6","IL8","IL10",
                                     "IL12A","IL12B","IL12RB1","IL12RB2","IL13","IL13RA1","IL13RA2","IL18","IFNA1","IFNB1",
                                     "IFNG","IRF1","IRF2"),]

############################
# extract barcodes for microarray and rna_seq 1098

rna_bars <- colnames(brca_data.rna_seq)
micr_bars <- colnames(brca_data.mrnaArray)

#rna_rownames <- row.names(brca_data.mrnaArray)
#rna_rownames <- row.names(brca_data.rna_seq)

############################
# extract barcodes for all 1212 samples

tcytokines <- t(brca_cytokines)  # transpose brca_cytokines
cytokines_barcode <- tcytokines["bcr_patient_barcode",]
cytokines_short_barcode <- tcytokines["barcode_short",]


colnames(tcytokines) <- cytokines_short_barcode[1:ncol(tcytokines)]
rna_cytokines <- tcytokines[-1,]
rna_cytokines_short_bcr <- tcytokines[-1,]
dim(rna_cytokines_short_bcr)
rna_cytokines_short_bcr <- rna_cytokines_short_bcr[-nrow(rna_cytokines_short_bcr),]  #trims last(tail) row

###############################
# match er+/er- barcodes and extract from my_data

col.num <- which(colnames(rna_cytokines_short_bcr) %in% er_pos_bcr2)
rna_cytokines_pos <- rna_cytokines_short_bcr[,col.num]
dim(rna_cytokines_pos)

er_pos_bcr2 <- colnames(rna_cytokines_pos)
#er_pos_cytokines <- rna_cytokines_short_bcr[,sort(c(col.num, col.num - 1))]
#dim(er_pos_cytokines)

col.num2 <- which(colnames(rna_cytokines_short_bcr) %in% er_neg_bcr2)
rna_cytokines_neg <- rna_cytokines_short_bcr[,col.num2]
dim(rna_cytokines_neg)

er_neg_bcr2 <- colnames(rna_cytokines_neg)
#er_neg_cytokines <- rna_cytokines_short_bcr[,sort(c(col.num2, col.num2 - 1))]
#dim(er_neg_cytokines)


###################################
# get sample types barcodes and matched types from TCGA


# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
#SS <- TCGAquery_SampleTypes(rna_bars,c("TP","NT"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(rna_bars,c("NT","TP"))
rna2_matched <- TCGAquery_MatchedCoupledSampleTypes(rna2,c("NT","TP"))

# retrieve matched samples in rna and microarray barcode groups
rna_matched <- TCGAquery_MatchedCoupledSampleTypes(rna_bars,c("NT","TP"))
micr_matched <- TCGAquery_MatchedCoupledSampleTypes(micr_bars,c("NT","TP"))

# check for matching barcodes in rna_matched and micr_matched
rna_micr_matched<- c()    # create a vector of all matching barcodes across rna and mircoarray vectors
for (i in rna_matched) {
  for (j in micr_matched) {
    if(i == j)
      rna_micr_matched<- append(rna_micr_matched,i) 
  }
} 

# rna matched normal and tumor
rna_matched_tumor <- TCGAquery_SampleTypes(rna_matched,"TP")   # primary tumor
rna_matched_normal <- TCGAquery_SampleTypes(rna_matched,"NT")  # normal tissue
rna2_normal <- TCGAquery_SampleTypes(cytokines_long_bcr,"NT")  # normal tissue bigger dataset
rna_matched_tumor<- sort(rna_matched_tumor)
rna_matched_normal<- sort(rna_matched_normal)

rna_matched_TN_TP_df <- data.frame(rna_matched_normal,rna_matched_tumor)

# microarray matched normal and tumor
micr_matched_tumor <- TCGAquery_SampleTypes(micr_matched,"TP")
micr_matched_normal <- TCGAquery_SampleTypes(micr_matched,"NT")
micr_matched_tumor<- sort(micr_matched_tumor)
micr_matched_normal<- sort(micr_matched_normal)

micr_matched_TN_TP_df <- data.frame(micr_matched_normal,micr_matched_tumor)

SSS_sorted <- sort(SSS)
SSS_sorted

#############################

rna_cytokines_pos_neg <- ER_classify(rna_bars)
micr_cytokines_pos_neg <- ER_classify(micr_bars)

er_pos_df <- data.frame(er_pos_cytokines)
er_neg_df <- data.frame(er_neg_cytokines)

col.num3 <- which(colnames(rna_cytokines_) %in% rna_micr_matched)
matched_rna_marray_cytokines <- rna_cytokines[,sort(c(col.num3, col.num3 - 1))]


############################
# Analysis

#find tumor normal in microarray
#remove normal from rna_seq
#remove normal from microarray

#rna_Seq analysis- deseq2

#- rna_seq data analysis er pos / er neg

#- microarray data (ER_classify to er+/er-)

#rna_seq
#- matched normal vs matched tumor er +ve
#- matched normal vs matched tumor er -ve

#microArray
#- matched normal vs matched tumor er +ve
#- matched normal vs matched tumor er -ve

# remove normal from cytokines_long_bcr  1212-112=1,100 and er_pos & er_neg 1149
rna2_cytokines_tumor_only <- cytokines_long_bcr[!cytokines_long_bcr %in% rna2_normal]
er_pos_tumor_only <- rna2_normal[!er_pos_bcr2 %in% rna2_normal]
er_neg_tumor_only <- rna2_normal[!er_neg_bcr2 %in% rna2_normal]

#- larger rna_seq_cytokines set minus normal 
# (rna er+ve vs er-ve)
er_pos_tumor_bcr <- rna2_cytokines_tumor_only[!rna2_cytokines_tumor_only %in% er_neg_bcr2]
er_neg_tumor_bcr <- rna2_cytokines_tumor_only[!rna2_cytokines_tumor_only %in% er_pos_bcr2]

#- STRING algorithm on top regulated
# done
