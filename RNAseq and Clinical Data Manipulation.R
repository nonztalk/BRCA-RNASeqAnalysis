# This data is too large to upload on github
load("./data/RNAseq_and_Clinical.RData")

# -------------------- Deal with clinical data -------------------------
clinic_subtype <- clinic[, c("bcr_patient_barcode", "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc")]
clinic_subtype[] <- lapply(clinic_subtype, function(x) ifelse(grepl("\\[.*\\]", x), "", x))

# Operation on missing value
NAcount <- function(x) {
  return(sum(x == ""))
}
clinic_subtype$NAcount <- apply(clinic_subtype, 1, NAcount)
# Delete all rows have missing value larger than 2
# clinic_subtype <- subset(clinic_subtype, NAcount < 2)
# You can see there're also Equivocal and Indeterminate

# This is the most strict filtering condition
# Let's try delete all NA value
clinic_subtype <- subset(clinic_subtype, NAcount == 0)
# And delete Equivocal and Indeterminate records
clinic_subtype <- subset(clinic_subtype, her2_status_by_ihc != "Equivocal" & her2_status_by_ihc != "Indeterminate")
clinic_subtype <- subset(clinic_subtype, pr_status_by_ihc != "Indeterminate" & er_status_by_ihc != "Indeterminate")
# Generate subtype using combination of three signatures
subtypeName <- function(x) {
  return(paste0(unname(sapply(x, substring, 1, 1)), collapse = ""))
}
clinic_subtype$subType <- apply(clinic_subtype[, 2:4], 1, subtypeName)
# delete NPN and NPP
clinic_subtype <- subset(clinic_subtype, subType != "NPN" & subType != "NPP")

# --------------------- Deal with RNA-seq data ------------------------
# Choose only the patients have records in clinical data
col_with_clinic <- colnames(data_tumor)[-1][vapply(colnames(data_tumor)[-1], 
                          function(x) {
                            ifelse(substring(x, 1, 12) %in% clinic_subtype$bcr_patient_barcode, TRUE, FALSE)
                          }, logical(1))]
data_tumor_with_clinic <- data_tumor[, col_with_clinic]
row.names(data_tumor_with_clinic) <- data_tumor$gene_id
# Construct data frame for patient id in data_tumor_with_clinic and their 
# corresponding subtype
data_id_situation <- data.frame(patientID = colnames(data_tumor_with_clinic),
                                subtype = vapply(colnames(data_tumor_with_clinic),
                                                 function(x) {
                                                   return(clinic_subtype$subType[clinic_subtype$bcr_patient_barcode == substring(x, 1, 12)])
                                                 }, character(1)), stringsAsFactors = F)
# two subtype have only 8 and 2 records
table(data_id_situation$subtype)
# NNN NNP NPN NPP PNN PNP PPN PPP 
# 111  32   8   2  64  23 355  99
# I try to combine them together in DE
# data_id_situation$subtype[data_id_situation$subtype == "NPN"] <- "NPE"
# data_id_situation$subtype[data_id_situation$subtype == "NPP"] <- "NPE"

# ----------------------------- limma analysis ----------------------
library(limma)
# delete rows with all 0 records
data_tumor_with_clinic <- data_tumor_with_clinic[rowSums(data_tumor_with_clinic) > 0, ]
# log scale log(x+1)
data_tumor_with_clinic <- log(data_tumor_with_clinic + 1)

# limma process
SampleCategoryFactor <- factor(data_id_situation[, "subtype"])
design <- model.matrix(~-1+SampleCategoryFactor)
contrast <- apply(t(combn(colnames(design), 2)), 1, paste, collapse = "-")
contrast.matrix <- makeContrasts(contrasts = contrast, levels = design)
fit <- lmFit(data_tumor_with_clinic, design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit1)
dif <- topTable(fit2, coef = contrast, n = nrow(fit2), lfc = log2(2.5))
dif <- dif[dif$adj.P.Val < 0.001, ] # 1562 genes

# the final matrix
matrixModel <- t(data_tumor_with_clinic[row.names(data_tumor_with_clinic) %in% row.names(dif), ])
matrixModel <- cbind(matrixModel, 
      clinic_subtype[clinic_subtype$bcr_patient_barcode %in% substring(row.names(matrixModel), 1, 12), 2:4])
matrixModel$subtype <- apply(matrixModel[, 1358:1360], 1, subtypeName)






# heatmap
# do not show NPE
library(pheatmap)
data_id_situation <- data.frame(patientID = colnames(data_tumor_with_clinic),
                                subtype = vapply(colnames(data_tumor_with_clinic),
                                                 function(x) {
                                                   return(clinic_subtype$subType[clinic_subtype$bcr_patient_barcode == substring(x, 1, 12)])
                                                 }, character(1)), stringsAsFactors = F)

data_id_situation <- data_id_situation[order(data_id_situation$subtype), ]
plotMatrix <- as.data.frame(t(matrixModel[, 1:1357]))
plotMatrix <- plotMatrix[data_id_situation$patientID]

patientColor <- vapply(colnames(plotMatrix), function(x) {
  Colors <- RColorBrewer::brewer.pal(8, "Pastel1")
  if (data_id_situation[x, ]$subtype == "NNN") return(Colors[1])
  if (data_id_situation[x, ]$subtype == "NNP") return(Colors[2])
  if (data_id_situation[x, ]$subtype == "NPN") return(Colors[3])
  if (data_id_situation[x, ]$subtype == "NPP") return(Colors[4])
  if (data_id_situation[x, ]$subtype == "PNN") return(Colors[5])
  if (data_id_situation[x, ]$subtype == "PNP") return(Colors[6])
  if (data_id_situation[x, ]$subtype == "PPN") return(Colors[7])
  if (data_id_situation[x, ]$subtype == "PPP") return(Colors[8])
}, character(1))

plotMatrix <- as.matrix(plotMatrix)
annotation_col <- data.frame(subtype = factor(data_id_situation$subtype))
row.names(annotation_col) <- colnames(plotMatrix)
annotation_color <- list(setNames(nm = unique(data_id_situation$subtype), object = RColorBrewer::brewer.pal(8, "Pastel1")))
pheatmap::pheatmap(plotMatrix, color = colorRampPalette(c('blue', 'white','red'))(20),
                   cluster_cols = F, legend = T, annotation_col = annotation_col, annotation_legend = T,
                   annotation_colors = annotation_color, annotation_names_col = T, cluster_rows = T, show_colnames = F,
                   fontsize_row = 1)

