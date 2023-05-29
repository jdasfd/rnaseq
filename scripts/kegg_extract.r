#!/usr/bin/env Rscript

args <- commandArgs(T)

if (!grepl("hsa", args[1])) {
    stop("Input error: please input hsa number: hsa...")
}

suppressMessages(library(KEGGREST))
suppressMessages(library(EnrichmentBrowser))
suppressMessages(library(clusterProfiler))

# init
out <- paste0(args[1], ".lst")

# download the pathway
hsapathway <- downloadPathways("hsa")

# retrive gene sets
hsa <- getGenesets(org = "hsa", db = "kegg", gene.id.type = "SYMBOL", cache = TRUE, return.type="list")

# write gmt to a tmp file
writeGMT(hsa, gmt.file = "kegg_hsa.tmp.gmt")
kegg <- clusterProfiler::read.gmt("kegg_hsa.tmp.gmt")

keggsub <- subset(kegg, grepl(args[1], kegg$term))
genelist <- keggsub$gene

# write out result and empty the tmp file
write.table(genelist, file = out, quote = FALSE, row.names = FALSE, col.names = FALSE)
file.remove("kegg_hsa.tmp.gmt")
