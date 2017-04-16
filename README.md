# brca_TCGA_project
Effect of Cytokines on differential gene expression in ER pos vs ER neg breast cancer

## Introduction

Breast cancer is the second leading cause of cancer deaths among women. Hereditary and environmental factors are known to predispose women to breast cancer. Genes and cytokines in the body are also known to play a role in the development or defense in the body against breast cancer. The growth of the tumor responds to stimulation or inhibition by these factors. We sought to investigate the differential expression of certain cytokines across estrogen receptor positive and estrogen receptor negative breast cancer subtypes and between normal and tumor tissue samples in women with breast cancer. Microarray and RNASeq analyses were used to quantify gene expression levels for genomic data on 1,098 tissue samples obtained from The Cancer Genome Atlas (TCGA) Data Portal (now absorbed into the Genomic Data Commons Project).  From these analyses, we are able to contribute to the current understanding of the biological factors that influence the development or treatment of breast cancer.           

Cytokines are any of a number of substances, such as interferon, interleukin, and growth factors, that are secreted by certain cells of the immune system and have an effect on other cells.

## Materials and Methods

The Cancer Genome Atlas (TCGA) is a collaboration between the National Cancer Institute (NCI) and the National Human Genome Research Institute (NHGRI) that has generated comprehensive, multi-dimensional maps of the key genomic changes in 33 types of cancer.  It is an open source data repository that holds genomic data (mutations, expression quantification, copy number variation) helps the cancer research community to improve the prevention, diagnosis, and treatment of cancer. In 2012, TCGA and other similar NIH funded open source projects were migrated to a single portal, the GDC data portal. 

The NCI's Genomic Data Commons (GDC) provides the cancer research community with a unified data repository that enables data sharing across cancer genomic studies in support of precision medicine.

A cancer is called estrogen-receptor-positive (or ER+) if it has receptors for estrogen. This suggests that the cancer cells, like normal breast cells, may receive signals from estrogen that could promote their growth. Testing for hormone receptors is important because the results help you and your doctor decide whether the cancer is likely to respond to hormonal therapy or other treatments. Hormonal therapy includes medications that either (1) lower the amount of estrogen in your body or (2) block estrogen from supporting the growth and function of breast cells. If the breast cancer cells have hormone receptors, then these medications could help to slow or even stop their growth. If the cancer is hormone-receptor-negative (no receptors are present), then hormonal therapy is unlikely to work.

Cytokines are any of a number of substances, such as interferon, interleukin, and growth factors, that are secreted by certain cells of the immune system and have an effect on other cells.

Using various tools (Firehose, TCGABiolinks, RTCGA Toolbox) we downloaded breast cancer datasets with clinical and genetic information from the TCGA/GDC portal into R.


## How we processed

Our analysis tool was R studio
Understanding the TCGA dataset (1098 rows, 1948 variables)
TCGA Barcodes, Patients, Samples, Aliquotes
Normal Tissue vs Primary Tumor
Estrogen receptor positive vs negative
Cytokines involved in Breast Cancer
Microarray vs RNASeq gene expression quantification data

## How we analyzed- Analysis tools in R
Limma
TCGA
Deseq
String

## How we got our results
Er + tumor vs normal
Er – tumor vs normal
Normal vs tumor
Er + vs Er -
Matched Er + tumor vs normal
Matched Er - tumor vs normal
Microarray vs RNAseq
Top Genes by logFC, p-value


