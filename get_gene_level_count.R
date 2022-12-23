#!/usr/bin/env Rscript
## This script quantifies gene level expression from Salmon/Kallisto

library(tximport)
args = commandArgs(trailingOnly=TRUE)

# specify the name of the quant.sf salmon file
quant_file = args[1]

# load table of transc \t gene correspondence
genes_file = args[2]
genes_dict = read.csv(genes_file, header=FALSE, sep='\t')

# provide type of analysis
analysis_type=args[3]

# get the output file from commandline
output = args[4]

#count_data = tximport(files=quant_file, type='salmon', tx2gene=genes_dict, ignoreTxVersion=TRUE)
#count_data = tximport(files=quant_file, type=analysis_type, tx2gene=genes_dict, ignoreTxVersion=TRUE)
#count_data = tximport(files=quant_file, type=analysis_type, tx2gene=genes_dict, countsFromAbundance="lengthScaledTPM")
count_data = tximport(files=quant_file, type=analysis_type, tx2gene=genes_dict)
write.table(count_data$counts, sep='\t', output, col.names=FALSE, quote=FALSE)

