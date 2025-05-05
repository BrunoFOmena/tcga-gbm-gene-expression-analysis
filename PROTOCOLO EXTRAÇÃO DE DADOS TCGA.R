
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "biomaRt"))


library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)


clinical_data <- read.delim("clinical.tsv", stringsAsFactors = FALSE)


barcodes <- unique(clinical_data$cases.submitter_id)


barcodes_30 <- barcodes[1:30]


query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = barcodes_30
)


GDCdownload(query)


data <- GDCprepare(query)

expr_counts <- assay(data)


expr_log2 <- log2(expr_counts + 1)


write.csv(expr_log2, file = "TCGA_GBM_counts_log2_normalized.csv")


genes_de_interesse <- c("EGFR", "TP53", "IDH1", "MGMT")


genes_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = genes_de_interesse,
                    mart = mart)

rownames(expr_log2) <- sub("\\..*", "", rownames(expr_log2))


expr_filtrado <- expr_log2[rownames(expr_log2) %in% genes_info$ensembl_gene_id, ]

write.csv(expr_filtrado, file = "genes_interesse_TCGA_GBM_log2.csv")
