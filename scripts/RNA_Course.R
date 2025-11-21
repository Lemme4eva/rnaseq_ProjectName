####For manual github setup ####
#renv::install("usethis)
#renv::install("gitcreds")
#renv::snapshot()
#library(usethis)
#library(gitcreds)
#gitcreds::gitcreds_set
#usethis::use_github

#Use other scripts
source("scripts/functions.R")
source("scripts/plotPCA2.R")
source("scripts/volcano.R")

#see information on a package
browseVignettes("DESeq2")

#install.packages("renv")
##renv::init()        # Initialize project environment
##renv::snapshot()    # Save package versions to renv.lock
##renv::restore()     # Recreate exact environment elsewhere


renv::install("here")
renv::install("fs")
renv::install("DT")
renv::install("ggplot2")
renv::install("org.Mm.eg.db")

renv::snapshot()

library(here)
library(fs)
library(ggplot2)

### Our practice data set

data_dir <- "." ## the current folder
matrix_file <- "GSE234563_Matrix.txt"

# download it ..
#download.file(
#  url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE234563&format=file&file=GSE234563%5FMatrix%2Etxt%2Egz",
#  destfile=matrix_file
#)


### .. read it:
M <- read.table(fs::path(data_dir, "data", "processed_data","counts_data",
                         matrix_file), #/Users/katharinalemberg/Documents/R/R_Course_MPI_advanced/rnaseq_ProjectName/data/processed_data
                header = TRUE,
                sep = "\t",
                quote = ""
)
# .. this is a huge matrix!
dim(M)

# .. display the first 6 rows (output next page):
dt <- DT::datatable(M |> head())
dt

#convert the column names to all-lowercase
colnames(M) <- tolower(colnames(M))

#discard genes with a total count < 10
min_count <- 10
M <- M[rowSums(M) >= min_count,]
#M <- M[(rowSums(M) >= min_count) &
#      apply(M,1,function(x) !any(x >= 10^3)) ,]

nrow(M)

#map matrix entries with fractional digits to the next higher integer value (also "round" ist possible)
M <- ceiling(M)

#### Setting Up the Input for DESeq2 ####
# Step 1: split the count matrix column names at "_"
#         (result: a list of vectors, where each vector holds the fragments extracted from one name)
colData <- strsplit(colnames(M),"_")

head(colData, n=2)

# Step 2: make a matrix, in which each vector from Step 1 is a row
colData <- Reduce(rbind,colData)

head(colData, n=4)

# Step 3: We know what the row and column names should be, so let's set them:
dimnames(colData) <- list(colnames(M),
                          c("age","sex","genotype","tissue","replicate")
)
head(colData, n=4)

# Step 4: The DESeq2 colData input should be a data.frame:
colData <- as.data.frame(colData)

# Step 5: Add a group which lists, for each sample, the levels of
#         all factors in this sample, concatenated by "." :
colData$group <-
  apply(colData,1,
        function(x) paste(x[c("age","sex","tissue","genotype")],
                          collapse="."))

options(width=150) ## increase output line length
# Show a bit more of the final result:
head(colData,n=4)



#### Is the nmrHas transgene effect strong enough to manifest irrespective of age, sex and tissue? ####
library(DESeq2)

dds_genotype <-
  # High level step 1: set up the object structure
  DESeqDataSetFromMatrix(countData = M,
                         colData = colData,
                         design= ~ 1 + genotype # in theory you could leave out the "1 + " because this is the default
  ) |>
  # High level step 2: "do the math" (--> Ali!)
  #                    size factors factors,
  #                    dispersion estimates,
  #                    statistical tests ...
  DESeq()

#available test results listed
resultsNames(dds_genotype)

#see genotype effect
res_genotype <-
  results(dds_genotype,
          name="genotype_nmrhas2_vs_creer"
  )


#put results in a nice table
dt <- DT::datatable(
  res_genotype |> as.data.frame() |>
    dplyr::filter(!is.na(padj), padj <= 0.05) |>
    dplyr::arrange(padj),

  options=list(scrollY="500px")
)

# convert numeric columns to scientific notation
dt <- dt |>
  DT::formatSignif(
    sapply(dt$x$data, is.numeric),
    digits=2
  )


#### Preparation for PCA ####
# The "blind" parameter tells the function NOT to use the input object's
# model formula for identifying "valid" (expected) variation.
# We do this here because we have not yet found an optimal formula.
# NOTE that if we have a trusted formula and want to use vst() to
# prepare an object for downstream analysis, then blind=FALSE should
# be used.
vsd_blind <- vst(dds_genotype,
                 blind=TRUE
)
saveRDS(vsd_blind,file="vsd_blind.rds")

plotPCA2(vsd_blind,
         intgroup = "group",
         ntop=500,
         PCs=c(x=1,y=2))

plotPCA2(vsd_blind,
         intgroup= "group",
         color=tissue,
         ntop=500,
         PCs=c(x=1,y=2), )


#Heatmap of the 500 Most Highly Expressed Genes
library(pheatmap)

selected <-
  order(rowMeans(counts(dds_genotype,
                        normalized=TRUE)
  ),
  decreasing=TRUE)[1:500]

df <- (colData(dds_genotype) |>
         as.data.frame())[,c("age",
                             "sex",
                             "genotype",
                             "tissue")]

pheatmap(assay(vsd_blind)[selected,],
         cluster_rows=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE,
         cluster_cols=TRUE,
         annotation_col=df)


#Boxplots of the Per-Sample Count Distributions : to see if the normalization "worked" - all medians should be aligned, but this is not the case, due to different distributions per tissue
## Reasoning: sometimes you do not know if you can analyze different parts separately
df <-
  rbind(reshape2::melt(log10(counts(dds_genotype,normalized=FALSE)+1)) |>
          mutate(mode="raw"),

        reshape2::melt(log10(counts(dds_genotype,normalized=TRUE)+1)) |>
          mutate(mode="normalized")
  ) |>
  # rename to make variable names more meaningful:
  dplyr::rename(gene=Var1, sample=Var2, count=value) |>

  # remove the trailing sample ID to create a common ID for all replicates
  # of a given factor level combination:
  mutate(group=sub("_\\d+$","",sample)) |>

  # ggplot needs to know the order of sub-plots ("raw" should come first)
  mutate(mode = factor(mode,
                       levels=c("raw","normalized")
  )
  )

ggplot(df,
       aes(y=sample,x=count,fill=group)
) +
  geom_boxplot(show.legend=FALSE) +
  facet_wrap(vars(mode),ncol=2) +
  theme_classic() +
  theme(axis.text=element_text(size=2))

#Boxplots of the Per-Sample Count Distributions
cnt <- counts(dds_genotype,normalized=TRUE)
medians_log10 <- log10(colMedians(cnt |> as.matrix())+1)

tapply(medians_log10, # what to group
       sub("_\\d+$","",
           names(medians_log10)
       ), # by what to group: sample name  without replicate ID
       mean  # function to apply to grouped values
) |>
  sort(decreasing=TRUE) # sort the result


#### DESeq for more factors, may expose more of a significant, although small, effects of the genotype factor ####
#shows much more results! Putting more factors in lets the system structure the data better, accounting for different variability
dds_additive <-
  # High level step 1: set up the object structure
  DESeqDataSetFromMatrix(countData = M,
                         colData = colData,
                         design= ~ age + sex + tissue + genotype
  ) |>
  DESeq()

resultsNames(dds_additive)


#Have a Look at the Genotype Effect in this Model
res_additive <-
  results(dds_additive, name="genotype_nmrhas2_vs_creer")

#make a table
dt <- DT::datatable(
  res_additive |> as.data.frame() |>
    dplyr::filter(!is.na(padj), padj <= 0.05) |>
    dplyr::arrange(padj),

  options=list(scrollY="50000px",
               lengthMenu = list(c(10, -1),
                                 c('10','All'))
  )
)

# convert numeric columns to scientific notation
dt <- dt |>
  DT::formatSignif(sapply(dt$x$data, is.numeric),
                   digits=2
  )
dt


###Group Contrasts
#new factors in colData which describe co-occurring levels of different factors in the same sample

options(width=150) ## increase output line length

## collapse (but don't remove!) columns age, sex, and tissue
## a new column "background":
colData$background <-
  apply(colData,1,
        function(x) paste(x[c("age","sex","tissue")],collapse="."))

head(colData,n=1)

#we express our desired pairwise comparisons in a table
contrasts <-
  data.frame(group= "group",
             A    = paste(colData[colData$genotype=="nmrhas2", "background"],
                          colData[colData$genotype=="nmrhas2", "genotype"],
                          sep="."
             ),
             B=     paste(colData[colData$genotype=="creer", "background"],
                          colData[colData$genotype=="creer", "genotype"],
                          sep="."
             )
  ) |> unique() # there is initially one copy per replicate!
saveRDS(contrasts,file="qmd_contrasts.rds")
head(contrasts,n=6)

#now run DESeq() with a 0 + group model. mean read count at a given “group” level is significantly different from zero
# then results() function, we can use its "contrast" argument to specifically compare a pair of “group” levels
#  it will inherit information on the full distribution from the DESeq() step.
# it takes longer to compute, but is a "less complex" model with less power because it only uses one factor, with 42 levels

dds_group <-
  DESeqDataSetFromMatrix(countData = M,
                         colData = colData,
                         design= ~ 0 + group) |>
  DESeq()

# a matrix allows to directly extract rows as vectors
# (unlike a data.frame):
m <- as.matrix(contrasts)

res_contrasts <- list()
for(i in 1:nrow(m)) {

  res_contrasts[[i]] <- results(dds_group,
                                contrast=m[i,]
  )
}
# the second name in each contrast pair without the trailing ".creer"
# is the background -- use as name:
names(res_contrasts)  <- sub("\\.creer$","",contrasts[,"B"])

res_contrasts[["old.female.intestine"]] |> as.data.frame() |>
  filter(padj <= 0.05) |> head(n=10)


#### Gene enrichment analysis ####
# Strategies to Manage Redundancy: Use tools that:
#  Summarize GO terms (e.g., REVIGO), Prune the DAG (e.g., parent-child approaches), Cluster enriched terms based on similarity
# Focus on most specific significant terms for clearer insights: packages like topGO, ClueGO, GO-Bayes, GO-Elit etc. help address hierarchy issues.

#using clusterProfiler

### Example of an ORA (Overrepresentation analysis)
#----------------------------------------------------------------------
# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(qs)

set.seed(420)

dds <- qread("dds_all_genotype.qs")
uni_genes <- results(dds) |> as.data.frame() |> # Universe is defined as genes that are actually expressed in the tissue
  filter(!is.na(padj)) |>
  rownames()
de_genes_up <- results(dds) |> as.data.frame() |>
  filter(!is.na(padj),padj <= 0.05, log2FoldChange > 0)|>
  rownames()
de_genes_dn <- results(dds) |> as.data.frame() |>
  filter(!is.na(padj),padj <= 0.05, log2FoldChange < 0)|>
  rownames()

ego_result_up <- enrichGO(gene = de_genes_up,
                          universe = uni_genes,
                          OrgDb = org.Mm.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)
# View results
#DT::datatable(head(ego_result_up))

# Visualize results
barplot(ego_result_up) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Barplot")

#----------------------------------------


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


#### Mini Project ####

###Exercise 1: Subset your data according to the tissue type assigned###

mat_spleen <- mat[, grepl("spleen", colnames(mat), ignore.case = TRUE)]

#Check data
dim(mat_spleen)       # rows x columns of the subset
colnames(mat_spleen)  # only spleen columns

####Exercise 2: Filter the samples to remove the low count genes.####
mat_spleen <- mat_spleen[rowSums(mat_spleen) >= 10, ]

####Exercise 3: Run DESeq2 package to get DE genes.####
coldata_spleen <- colnames(mat_spleen) |>
  strsplit("_") |>
  unlist() |>
  matrix(ncol = 5, byrow = TRUE) |>
  as.data.frame() |>
  setNames(c("age", "sex", "genotype", "tissue", "replicate"))

dds_spleen <- DESeqDataSetFromMatrix(
  countData = mat_spleen,
  colData = coldata_spleen,
  design = ~ age + sex + genotype
) |>
  DESeq()

qsave(dds_spleen, "results/differential_expression/dds_spleen.qs")
dds <- qread("results/differential_expression/dds_spleen.qs")

####Exercise 4: Use boxplot to see the distribution of the samples. What do you observe?####

df <-
  rbind(reshape2::melt(log10(counts(dds_spleen,normalized=FALSE)+1)) |>
          mutate(mode="raw"),

        reshape2::melt(log10(counts(dds_spleen,normalized=TRUE)+1)) |>
          mutate(mode="normalized")
  ) |>
  # rename to make variable names more meaningful:
  dplyr::rename(gene=Var1, sample=Var2, count=value) |>

  # remove the trailing sample ID to create a common ID for all replicates
  # of a given factor level combination:
  mutate(group=sub("_\\d+$","",sample)) |>

  # ggplot needs to know the order of sub-plots ("raw" should come first)
  mutate(mode = factor(mode,
                       levels=c("raw","normalized")
  )
  )

ggplot(df,
       aes(y=sample,x=count,fill=group)
) +
  geom_boxplot(show.legend=FALSE) +
  facet_wrap(vars(mode),ncol=2) +
  theme_classic() +
  theme(axis.text=element_text(size=2))

#Boxplots of the Per-Sample Count Distributions
cnt <- counts(dds_spleen,normalized=TRUE)
medians_log10 <- log10(colMedians(cnt |> as.matrix())+1)

tapply(medians_log10, # what to group
       sub("_\\d+$","",
           names(medians_log10)
       ), # by what to group: sample name  without replicate ID
       mean  # function to apply to grouped values
) |>
  sort(decreasing=TRUE) # sort the result

#### Exercise 5: Visualize the results using MA-plot and PCA.####

vsd_blind <- vst(dds_spleen,
                 blind=TRUE
)
saveRDS(vsd_blind,file="vsd_blind.rds")

plotPCA2(vsd_blind,
         intgroup = "genotype",
         ntop=500,
         PCs=c(x=1,y=2))

plotPCA2(vsd_blind,
         intgroup= "genotype",
         color=genotype,
         ntop=500,
         PCs=c(x=1,y=2))

plotPCA2(vsd_blind,
         intgroup= "sex",
         color=sex,
         ntop=500,
         PCs=c(x=1,y=2))

plotPCA2(vsd_blind,
         intgroup= "age",
         color=age,
         ntop=500,
         PCs=c(x=1,y=2))


## Use DESeq2's plotMA function
res_spleen <- results(dds_spleen, name="sex_male_vs_female")

plotMA(res_spleen, main="MA Plot Spleen", ylim=c(-10,10))


####Exercise 6: Save the results to a csv file####
res_spleen_df <- as.data.frame(res_spleen) |>
  rownames_to_column("gene")  # gene names in column

# Save to file
write.csv(res_spleen_df, file = "results/differential_expression/res_spleen_dds_results.csv", row.names = TRUE)


#### Exercise 7: Create a volcano plot to visualize the DE genes####

source("scripts/volcano.R")
p <- volcano(
  df = res_spleen_df,
  use_cols = c(2,4,5),  # specify columns
  lfc_thrsh = 1,   # e.g., log2 fold change > 1 or < -1
  p_thrsh   = 0.05,
  label.show = TRUE
)

p

##another volcano plot
# Extract results and convert to data frame

renv::install("tibble")
library(tibble)

res_spleen_df <- as.data.frame(res_spleen) |>
  rownames_to_column("gene")  # gene names in column

#add significance and threshold
lfc_thresh <- 1       # log2 fold change threshold
p_thresh   <- 0.05    # adjusted p-value threshold

#add significance column
res_spleen_df <- res_spleen_df |>
  mutate(
    significant = ifelse(!is.na(padj) & padj <= p_thresh & abs(log2FoldChange) >= lfc_thresh,
                         ifelse(log2FoldChange > 0, "Up", "Down"),
                         "NotSig")
  )

renv::install("ggrepel")
library(ggrepel)

# Select top 10 up and top 10 down
top_genes <- res_spleen_df |>
  filter(significant != "NotSig") |>
  arrange(padj) |>
  slice_head(n = 20)  # pick 20 smallest padj

ggplot(res_spleen_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_thresh), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10(p-value)", color = "Significance") +
  ggtitle("Volcano Plot with Top Genes")


#### Exercise 8: Perform GO enrichment (ORA analysis) on the DE genes(UP and Down). ####

uni_genes_spleen <- results(dds_spleen, name="sex_male_vs_female") |>
  as.data.frame() |>
  filter(!is.na(padj)) |>
  rownames()
de_genes_spleen_up <- results(dds_spleen, name="sex_male_vs_female") |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
  rownames()
de_genes_spleen_dn <- results(dds_spleen, name="sex_male_vs_female") |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange < 0) |>
  rownames()
de_genes_spleen_both <- results(dds_spleen, name="sex_male_vs_female") |>
  as.data.frame() |>
  filter(!is.na(padj), padj <= 0.05, log2FoldChange != 0) |>
  rownames()


ego_result_spleen_up <- enrichGO(
  gene = de_genes_spleen_up,
  universe = uni_genes_spleen,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

DT::datatable(head(ego_result_spleen_up)) # shows all genes that are now enriched for the pathways and how many there are

ego_result_spleen_dn <- enrichGO(
  gene = de_genes_spleen_dn,
  universe = uni_genes_spleen,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

DT::datatable(head(ego_result_spleen_dn))

# Visualizations
# x="Count"/"GeneRatio", color="p.adjust"/"qvalue", size="GeneRatio"/"Count"

dotplot(ego_result_spleen_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Dotplot")
dotplot(ego_result_spleen_dn, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Dotplot")

barplot(ego_result_spleen_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Barplot")
barplot(ego_result_spleen_dn, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Barplot")

cnetplot(ego_result_up, categorySize="pvalue", foldChange=NULL) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot")
cnetplot(ego_result_dn, categorySize="pvalue", foldChange=NULL) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Cnetplot")

ego_result_spleen_up <- enrichplot::pairwise_termsim(ego_result_spleen_up)
ego_result_spleen_dn <- enrichplot::pairwise_termsim(ego_result_spleen_dn)
emapplot(ego_result_spleen_up, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Upregulated Genes : Emapplot")
emapplot(ego_result_spleen_dn, showCategory=10) +
  ggtitle("GO Enrichment Analysis - Downregulated Genes : Emapplot")


#### Exercise 9: Perform GSEA and use dot plot to visualize ####


#Get DESeq2 results
res_spleen_sex <- results(dds_spleen, name = "sex_male_vs_female")

##Create a ranked gene list for GSEA
#use the Wald statistic (res$stat):

df <- res_spleen_sex |>
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








#Heatmap of the 500 Most Highly Expressed Genes
library(pheatmap)

selected <-
  order(rowMeans(counts(dds_genotype,
                        normalized=TRUE)
  ),
  decreasing=TRUE)[1:500]

df <- (colData(dds_genotype) |>
         as.data.frame())[,c("age",
                             "sex",
                             "genotype",
                             "tissue")]

pheatmap(assay(vsd_blind)[selected,],
         cluster_rows=FALSE,
         show_rownames=FALSE,
         show_colnames=FALSE,
         cluster_cols=TRUE,
         annotation_col=df)






