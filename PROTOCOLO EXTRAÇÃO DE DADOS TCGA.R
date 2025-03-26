# Instalar pacotes se ainda não tiver
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "biomaRt"))

# Carregar os pacotes
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)

# 1. Lê o arquivo clinical.tsv
clinical_data <- read.delim("clinical.tsv", stringsAsFactors = FALSE)

# 2. Extrai os barcodes (removendo duplicados)
barcodes <- unique(clinical_data$cases.submitter_id)

# 3. Seleciona 30 barcodes
barcodes_30 <- barcodes[1:30]

# 4. Consulta os dados de contagem (STAR - Counts)
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = barcodes_30
)

# 5. Faz o download dos dados
GDCdownload(query)

# 6. Prepara os dados para análise
data <- GDCprepare(query)

# 7. Extrai a matriz de contagens
expr_counts <- assay(data)

# 8. Normaliza com log2(count + 1)
expr_log2 <- log2(expr_counts + 1)

# 9. Salva a matriz completa (normalizada) como CSV
write.csv(expr_log2, file = "TCGA_GBM_counts_log2_normalized.csv")

# 10. Define genes de interesse
genes_de_interesse <- c("EGFR", "TP53", "IDH1", "MGMT")

# 11. Converte nomes para Ensembl IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = genes_de_interesse,
                    mart = mart)

# ✅ CORRIGE os Ensembl IDs com sufixo (ex: .14, .15)
rownames(expr_log2) <- sub("\\..*", "", rownames(expr_log2))

# 12. Filtra a matriz só com os genes desejados
expr_filtrado <- expr_log2[rownames(expr_log2) %in% genes_info$ensembl_gene_id, ]

# 13. Salva o CSV final só com os genes de interesse
write.csv(expr_filtrado, file = "genes_interesse_TCGA_GBM_log2.csv")
