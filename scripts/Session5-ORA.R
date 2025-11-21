#### Day 2 ####
# Load necessary libraries

renv::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(qs)
library(DT)
library(enrichplot)


#Use other scripts
source("scripts/functions_sessions2_4.R")
source("scripts/volcano.R")

#### ORA (overrepresentation analysis) ####
mat <- read.table(fs::path(data_dir, "data", "processed_data","counts_data",
                         matrix_file))
colnames(mat) <- tolower(colnames(mat))
mat <- ceiling(mat)
mat <- mat[rowSums(mat) >= 10, ]

coldata <- colnames(mat) |>
    strsplit("_") |>
    unlist() |>
    matrix(ncol = 5, byrow = TRUE) |>
    as.data.frame() |>
    setNames(c("age", "sex", "genotype", "tissue", "replicate"))

dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = coldata,
    design = ~ age + sex + tissue + genotype
) |>
    DESeq()

qsave(dds, "results/differential_expression/dds_all_genotype.qs")
dds <- qread("results/differential_expression/dds_all_genotype.qs")
uni_genes <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj)) |>
    rownames()
de_genes_up <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
    rownames()
de_genes_dn <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj), padj <= 0.05, log2FoldChange < 0) |>
    rownames()
de_genes_both <- results(dds) |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange != 0) |>
  rownames()


ego_result_up <- enrichGO(
    gene = de_genes_up,
    universe = uni_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

DT::datatable(head(ego_result_up)) # shows all genes that are now enriched for the pathways and how many there are

ego_result_dn <- enrichGO(
    gene = de_genes_dn,
    universe = uni_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

DT::datatable(head(ego_result_dn))

# Visualizations
# x="Count"/"GeneRatio", color="p.adjust"/"qvalue", size="GeneRatio"/"Count"

dotplot(ego_result_up, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Dotplot")
dotplot(ego_result_dn, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Dotplot")

barplot(ego_result_up, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Barplot")
barplot(ego_result_dn, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Barplot")

cnetplot(ego_result_up, categorySize="pvalue", foldChange=NULL) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot")
cnetplot(ego_result_dn, categorySize="pvalue", foldChange=NULL) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Cnetplot")

ego_result_up <- enrichplot::pairwise_termsim(ego_result_up)
ego_result_dn <- enrichplot::pairwise_termsim(ego_result_dn)
emapplot(ego_result_up, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Emapplot")
emapplot(ego_result_dn, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Emapplot")



## does not work for me, probably dplyr select is masked by DESeq2
#result_up <- results(dds) |>
#	as.data.frame() |>
#	filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
#	select(log2FoldChange) |>
#	tibble::rownames_to_column(var = "gene") |>
#	pull(log2FoldChange,gene)

#adding the p value -- then you can see which genes are more significant in the plot
result_up <- results(dds) |>
  as.data.frame() |>
  dplyr::filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
  dplyr::select(log2FoldChange) |>
  tibble::rownames_to_column(var = "gene") |>
  dplyr::pull(log2FoldChange, gene)

cnetplot(ego_result_up, foldChange = result_up) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot with log2FC")


#### Gene Set Enrichment Analysis ####
# GSEA gene lists are ranked, by a principally arbitrary rank metric, reflecting the interest of the experimenter
# the rank will be "walked down" and overlap with the gene set will increase the score, and the score will decrease again
# for non-overlapping genes, giving a "curve", thus also giving a "spacial" information, in which rank area the overlappings are

#using e.g. clusterProfiler::gseGO (uses fgsea algorithm)

#Get DESeq2 results
res <- results(dds)

##Create a ranked gene list for GSEA
#use the Wald statistic (res$stat):

df <- res |>
  as.data.frame() |>
  dplyr::filter(!is.na(padj)) |>
  dplyr::arrange(desc(stat))

geneList <- df$stat
names(geneList) <- rownames(df)

gsea_results <- runGSEA(
  geneList = geneList,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",   # or "ENSEMBL", depending on your data
  ontologies = c("BP","MF","CC")
)


#1. Check what’s in the list
names(gsea_results)

#2. Convert to a data frame
bp_df <- as.data.frame(gsea_results$BP)
mf_df <- as.data.frame(gsea_results$MF)
cc_df <- as.data.frame(gsea_results$CC)

#3. View top pathways
head(bp_df[order(bp_df$p.adjust), ])

#4. Extract only the significant pathways
sig_bp <- bp_df[bp_df$p.adjust <= 0.05, ]

#5. Access specific info
# The column with actual contributing genes is 'core_enrichment'
bp_df$core_enrichment[1]
# Split into vector
genes <- strsplit(bp_df$core_enrichment[1], "/")[[1]]


sig_bp[, c("Description", "NES", "p.adjust")] #GO term descriptions and NES

#6 Plot results
bp_df <- as.data.frame(gsea_results$BP)
bp_df <- bp_df[order(bp_df$p.adjust), ][1:10, ]  # top 10

ggplot(bp_df, aes(x = reorder(Description, NES), y = NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Normalized Enrichment Score (NES)", x = "GO Term")

# Dot plot
dotplot(gsea_results$BP, showCategory = 15)

# Enrichment plot for one gene set
gseaplot2(gsea_results$BP, geneSetID = sig_bp$ID[1])


