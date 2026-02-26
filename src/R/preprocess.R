# src/R/preprocess.R
# Title: Baron Pancreas Harmonization (Mouse to Human)
# Author: Gemini-Co-Developer

library(Seurat)
library(biomaRt)
library(dplyr)

cat("🧬 Starting Phase 2: Ortholog Mapping...\n")

# 1. Connect to Ensembl (The "Dictionary" for Gene Mapping)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 2. Define the Mapping Function
get_orthologs <- function(mouse_genes) {
  getLDS(attributes = c("mgi_symbol"), 
         filters = "mgi_symbol", 
         values = mouse_genes, 
         mart = mouse, 
         attributesL = c("hgnc_symbol"), 
         martL = human)
}

# 3. Simulated Baron Data Example 
# (In production, this loads from data/raw/baron_mouse.csv)
mouse_genes_example <- c("Ins1", "Ins2", "Gcg", "Sst", "Ppy")

cat("🔍 Mapping Mouse symbols to Human HGNC...\n")
mapping_table <- get_orthologs(mouse_genes_example)

# 4. Display Results
print(mapping_table)

# 5. Export for PyTorch Phase
# write.csv(mapping_table, "data/processed/gene_mapping.csv", row.names = FALSE)
cat("✅ Phase 2 Complete: Genes mapped and ready for alignment.\n")