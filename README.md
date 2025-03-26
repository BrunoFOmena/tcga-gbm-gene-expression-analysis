# Análise de Expressão Gênica em Glioblastoma Multiforme (TCGA-GBM)

Este repositório contém um pipeline de pré-processamento de dados de expressão gênica do projeto TCGA-GBM, com foco na normalização dos dados e extração de genes de interesse relacionados ao Glioblastoma Multiforme.

## Objetivo

O objetivo principal deste projeto é realizar o download, a normalização e a filtragem dos dados de expressão gênica de amostras tumorais de glioblastoma disponíveis no The Cancer Genome Atlas (TCGA), visando facilitar análises posteriores, como identificação de biomarcadores ou validação de alvos terapêuticos.

## Genes de Interesse

Os seguintes genes foram selecionados com base em sua relevância descrita na literatura sobre glioblastoma:

- **EGFR**
- **TP53**
- **IDH1**
- **MGMT**

## Etapas Realizadas

1. **Leitura dos dados clínicos**  
   Importação do arquivo `clinical.tsv` com as informações das amostras.

2. **Seleção de barcodes**  
   Seleção de 30 barcodes únicos para amostras do tipo tumor primário.

3. **Download dos dados de expressão**  
   Consulta e download dos dados de contagem de RNA-Seq (workflow: STAR - Counts) usando o pacote `TCGAbiolinks`.

4. **Normalização**  
   Aplicação de transformação log2(count + 1) para normalizar os dados.

5. **Mapeamento de genes de interesse**  
   Conversão dos símbolos gênicos para IDs Ensembl utilizando o pacote `biomaRt`.

6. **Filtragem da matriz**  
   Seleção da matriz apenas com os genes de interesse e exportação para CSV.

## Requisitos

- R >= 4.0
- Pacotes:
  - `TCGAbiolinks`
  - `SummarizedExperiment`
  - `biomaRt`

## Execução

Para executar o pipeline, basta rodar o script contido em `scripts/pre_processamento_TCGA_GBM.R`. O script realiza automaticamente todas as etapas descritas acima.
