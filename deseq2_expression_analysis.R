#!/usr/bin/env Rscript
#
## This script runs DESeq2 analysis for requested type between 2 requested species

## Install required packeges if needed
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("argparse")

## Load argument parser package
suppressPackageStartupMessages(library(argparse))

## Initiate argument parser
parser = ArgumentParser()
parser$add_argument("-w", "--workdir", type="character", help="abs path to work (current) dir")
parser$add_argument("-d", "--dircounts", type="character", help="relative path to directory with all count files by sample")
parser$add_argument("-c", "--counts", type="character", help="table (.tsv) with all counts for all samples")
parser$add_argument('-m', "--meta", type="character", help="meta data file: all info about samples")
parser$add_argument('-s1', '--sp1', type='character', help='TEST species for DE analysis')
parser$add_argument('-s2', '--sp2', type='character', help='REFERENCE species for DE analysis')
parser$add_argument('-t', '--type', type="character", help="type, e.g. liver / pectoralis / heart; OR liver,pectoralis if running multifactorial")
parser$add_argument('-od', '--outdeseq', type="character", help="deseq2 output file: gene, foldChange, padj..")
parser$add_argument('-on', '--outnorm', type="character", help="output file with norm counts for requested samples")
parser$add_argument('-multi', '--multi', action='store_true', default=FALSE, help='specify if want to run multifactorial analysis')
args = parser$parse_args()

## Load other required packages
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))

## Parse arguments from command line
curr_dir = args$workdir
setwd(curr_dir)
counts_dir = args$dircounts
counts_file = args$counts
meta_file = args$meta
sp_target = args$sp1
sp_ref = args$sp2
t = args$type
out_file_deseq = args$outdeseq
out_norm_counts = args$outnorm
multi = args$multi

## Load all counts
all_sample_dirs = list.dirs(path=counts_dir, full.names=FALSE, recursive=FALSE)
all_count_files = paste(counts_dir, all_sample_dirs, counts_file, sep='/')

for (i in 1:length(all_count_files)){
  f = all_count_files[i]
  l_name = all_sample_dirs[i]
  
  # if the merged dataset doesn't exist, create it
  if ( ! exists("df_merged") && file.exists(f)) {
    df_merged = read.csv(f, header=FALSE, sep="\t")
    colnames(df_merged) = c('gene', l_name)
  }
  
  # if the merged dataset does exist, append to it
  else if (exists("df_merged") && file.exists(f)) {
    temp_df = read.csv(f, header=FALSE, sep="\t")
    colnames(temp_df) = c('gene', l_name)
    df_merged = merge(df_merged, temp_df, by='gene', all=TRUE)
    rm(temp_df)
  }
}

## Round counts
counts_table_round <- df_merged %>% mutate_if(is.numeric, round)

## Write ALL count ALL gene table into a file
# all_lib_all_gene_file = 'ALL_lib.ALL_gene_table.tsv'
# write.table(counts_table_round, sep="\t", all_lib_all_gene_file, quote=FALSE)
# quit(status=0)

## Load sample info
meta_data = read.csv(meta_file, header=TRUE, sep='\t')

## Set up DESeq2 analysis
run_deseq_analysis = function(meta, counts, sources, sp_ref, multi) {
  ## Get data for the requested type
  meta_sources = meta[meta$lib %in% sources, ]
  counts_sources = na.omit(counts[sources])
  
  ## Run DESeq2 analysis
  if (multi){
    ## if multifactorial
    dds = DESeqDataSetFromMatrix(countData=counts_sources, 
                                 colData=meta_sources, 
                                 design=~species+species:type,
                                 tidy=TRUE)
  } else {
    ## if pairwise
    dds = DESeqDataSetFromMatrix(countData=counts_sources, colData=meta_sources, design=~species, tidy=TRUE)
  }
  
  # set the factor level
  dds$species = relevel(dds$species, ref = sp_ref)
  dds = DESeq(dds)
  return(dds)
}


### Run DESeq2 analysis

## Set up
if (multi){
  ## if multifactorial
  types = strsplit(t, split=',', fixed=TRUE)[[1]]
  t_ref = types[1]
  t_target = types[2]
  print('Running multifactorial analysis . . .')
  lib_type = subset(meta_data, (meta_data$type == t_target | meta_data$type == t_ref) & 
                      (meta_data$species == sp_target | meta_data$species == sp_ref))$lib
} else {
  ## if pairwise
  print('Running pairwise analysis . . .')
  lib_type = subset(meta_data, meta_data$type == t & (meta_data$species == sp_target | meta_data$species == sp_ref))$lib
  
}

sources_type = c('gene', as.vector(lib_type)) 
dds = run_deseq_analysis(meta_data, counts_table_round, sources_type, sp_ref, multi)


### Process deseq2 output and write to output files

# Write normalized counts to a file
norm_counts_res = counts(dds, normalized=TRUE)
write.table(norm_counts_res, sep="\t", out_norm_counts, quote=FALSE)

if (multi){
  ## Multifactor results
  mod_mat <- model.matrix(design(dds), colData(dds))
  # Define coefficient vectors for each condition
  test_test <- colMeans(mod_mat[dds$species == sp_target & dds$type == t_target, ])
  test_ref  <- colMeans(mod_mat[dds$species == sp_target & dds$type == t_ref, ])
  ref_test  <- colMeans(mod_mat[dds$species == sp_ref    & dds$type == t_target, ])
  ref_ref   <- colMeans(mod_mat[dds$species == sp_ref    & dds$type == t_ref, ])
  
  res_ref <- results(dds, contrast = ref_test - ref_ref)
  res_interaction <- results(dds, contrast = (test_test - test_ref) - (ref_test - ref_ref))
  
  # Order by adjusted p-value
  res_ref_oredered = res_ref[order(res_ref$padj), ]
  res_interaction_oredered = res_interaction[order(res_interaction$padj), ]
  # Write results to files
  out_file_ref = paste0('ref.', out_file_deseq)
  out_file_interaction = paste0('interaction.', out_file_deseq)
  write.table(res_ref_oredered, sep="\t", out_file_ref, quote=FALSE)
  write.table(res_interaction_oredered, sep="\t", out_file_interaction, quote=FALSE)
  
} else {
  ## Pairwise results
  res = results(dds)
  res_ordered = res[order(res$padj), ]
  write.table(res_ordered, sep="\t", out_file_deseq, quote=FALSE)
}



