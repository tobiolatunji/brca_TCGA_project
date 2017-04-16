#' classify the patients into ER+ and ER-
#'
#' @param barcodes A vector of barcodes which can either be a 12-bit version or long version
#' @return A list consists of the clinical.patient data frame 
#' alone with two vectors of barcodes stand for ER+ and ER- respective
#' @examples 
#' barcodes <-c("TCGA-XX-A89A", "TCGA-Z7-A8R5", "TCGA-Z7-A8R6")
#' ER_classify(barcodes)
ER_classify <- function(barcodes){
  if(!is.null(barcodes) && length(barcodes)>0){
    if( nchar(barcodes[1])>12){
      barcodes<-substr(barcodes,1,12)
    }
    clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical",barcode = barcodes)
    tryCatch(GDCdownload(clin.query), error = function(e) GDCdownload(clin.query, method = "client"))
    clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
    er_pos<-TCGAquery_clinicFilt( barcodes, clinical.patient, ER="Positive")
    er_neg<-TCGAquery_clinicFilt( barcodes, clinical.patient, ER="Negative")
    return ( list(clinical.patient, er_pos,er_neg))
  }
}


