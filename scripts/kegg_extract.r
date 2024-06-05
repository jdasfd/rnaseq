#!/usr/bin/env Rscript

# Env dependence
# R (4.2.3)
# Bioconductor (3.18)

# Usage
# Rscript kegg_extract.r <hsa_number>

# If there is anything wrong, please following the steps here:

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.19")
# 
# packages <- c("KEGGREST", "EnrichmentBrowser", "clusterProfiler")
# BiocManager::install(setdiff(packages, rownames(installed.packages())))

args <- commandArgs(T)

if (!grepl("hsa", args[1])) {
    stop("Input error: please input hsa number: hsa...")
}

suppressMessages(library(KEGGREST))
suppressMessages(library(EnrichmentBrowser))
suppressMessages(library(clusterProfiler))

# init
out <- paste0(args[1], ".txt")

# download the pathway
hsapathway <- downloadPathways("hsa")

# retrive gene sets
hsa <- getGenesets(org = "hsa", db = "kegg", gene.id.type = "SYMBOL", cache = TRUE, return.type="list")

# write gmt to a tmp file
writeGMT(hsa, gmt.file = "kegg_hsa.tmp.gmt")
kegg <- clusterProfiler::read.gmt("kegg_hsa.tmp.gmt")

# extract the subset according to the hsa accession
# if you want to directly search this, please replace args[1] to "hsa+numer"
keggsub <- subset(kegg, grepl(args[1], kegg$term))
genelist <- keggsub$gene

# write out result and empty the tmp file
# if you do not want this step, skip write table or just use your own way to export
write.table(genelist, file = out, quote = FALSE, row.names = FALSE, col.names = FALSE)
file.remove("kegg_hsa.tmp.gmt")
