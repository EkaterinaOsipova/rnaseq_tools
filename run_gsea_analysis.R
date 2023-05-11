#!/usr/bin/env Rscript

## GSEA ---------------------------------------------------------

## Load argument parser package
suppressPackageStartupMessages(library(argparse))

## Initiate argument parser
parser = ArgumentParser()
parser$add_argument("-w", "--workdir", type="character", help="abs path to work (current) dir")
parser$add_argument("-g", "--genelist", type="character", help="table (.tsv) gene - effectSize")
parser$add_argument("-o", "--outfile", type="character", help="GSEA results output file name (raw p-value cutoff <.01)")
args = parser$parse_args()

## Load clusterProfiler and related packages
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))
# install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(tibble)

## Parse arguments from command line
curr_dir = args$workdir
setwd(curr_dir)
file_name = args$genelist
gsea_out_file = args$outfile


## Load gene symbols and weights
file_df <- read.csv(file_name, header=FALSE, sep='\t')
colnames(file_df) <- c('gene', 'effect')
gene_dos <- na.omit(file_df[c('gene', 'effect')])
gene_dos <- gene_dos[order(-gene_dos$effect), ]
gene_dos <- gene_dos[gene_dos$effect != 0, ]


## Perform GSEA with clusterProfiler
gene_list <- gene_dos$effect
names(gene_list) = gene_dos$gene

gse <- gseGO(geneList         = gene_list,
             OrgDb         = org.Hs.eg.db,
             keyType       = "SYMBOL",
             ont           = "BP",
             pAdjustMethod = 'BH',
             pvalueCutoff  = 0.05,
             minGSSize = 40,
             maxGSSize = 200)


gse_result = gse@result
write.table(gse_result[gse_result$pvalue < 0.01, ], row.names=F, quote=F, file=gsea_out_file, sep="\t")

