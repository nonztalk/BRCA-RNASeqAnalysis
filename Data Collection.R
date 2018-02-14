source("./TCGA-Assembler/Module_A.R")
source("./TCGA-Assembler/Module_B.R") # if just download the data, Module B is not necessary

# download BRCA RNA-Seq data
DownloadRNASeqData(
  cancerType = "BRCA", # set cancer type
  assayPlatform = "gene.normalized_RNAseq", # Normalized gene expression quantification
  tissueType = "TP", # Primary Solid Tumor
  saveFolderName = "./data" # where to save this file
)

# download clinical data
DownloadBiospecimenClinicalData(
  cancerType = "BRCA",
  saveFolderName = "./data"
)