library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Your list of genes
genes_to_test <- proof_genes
genes_to_test <- toupper(trimws(genes_to_test))

# Convert Symbols to Entrez IDs
gene_conv <- bitr(genes_to_test, 
                  fromType = "SYMBOL", 
                  toType   = "ENTREZID", 
                  OrgDb    = org.Hs.eg.db)

# Run enrichment using the Entrez IDs
go_results <- enrichGO(gene = gene_conv$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = 'ENTREZID', 
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 1,
                       readable = TRUE) #
# 3. View the results
head(go_results, n = 10)


  as.data.frame(ego_mf_sim)

