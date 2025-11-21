library(dplyr)
library(ggplot2)
library(pheatmap)
library(DESeq2)
source("R/functions.R")
source("R/plotPCA2.R")

data_dir <- "." ## the current folder
matrix_file <- "GSE234563_Matrix.txt.gz"

M <- read.table(fs::path(data_dir,
                         matrix_file),
                header = TRUE,
                sep = "\t",
                quote = ""
                )
dim(M)
dt <- DT::datatable(M |> head(),options(scrollX=TRUE))
dt

##First of all, we **convert the column names to all-lowercase**, to avoid case issues:
colnames(M) <- tolower(colnames(M))

##Then we **discard genes with a total count < 10**, as was done in the paper: 

min_count <- 10
M <- M[rowSums(M) >= min_count,]

nrow(M)
##Finally we **map matrix entries with fractional digits to the next higher integer value**\        
##(the differential expression analysis expects integer counts):

M <- ceiling(M)

## Introducing the DESeq2 Package
library(DESeq2)

browseVignettes("DESeq2")

## Setting Up the Input for DESeq2


# We do not have the colData -- but it is actually easy to create this matrix starting from the count matrix column names:


# Step 1: split the count matrix column names at "_"
#         (result: a list of vectors, where each vector holds the fragments extracted from one name)
colData <- strsplit(colnames(M),"_") 
head(colData, n=2) 

colData <- Reduce(rbind,colData)
head(colData, n=2) 

## Setting Up the Input for DESeq2

# Step 3: We know what the row and column names should be, so let's set them:
dimnames(colData) <- list(colnames(M),
c("age","sex","genotype","tissue","replicate")
                          )
head(colData, n=2) 

# Step 4: The DESeq2 colData input should be a data.frame:
colData <- as.data.frame(colData)


##colData$group <- 
##  purrr::reduce(colData[,c("age","sex","genotype","tissue")],
##                paste, sep=".")
apply(colData[,c("age","sex","genotype","tissue")],1,paste,collapse=".")

# Show a bit more of the final result:
head(colData, n=9) 

## A Bold Try: Model the Data With Genotype Only

##The most simple question we can ask to our data matrix is: \           
##**Is the nmrHas transgene effect strong enough to manifest irrespective of age, sex and tissue?**  
##Starting from this question, we will dig deeper both into the dataset and into modeling.


library(DESeq2)

dds_genotype <-
  # High level step 1: set up the object structure
    DESeqDataSetFromMatrix(countData = M,
                           colData = colData,
                           design= ~ 1 + genotype 
                           ) |>
  # High level step 2: "do the math" (--> Ali!)
  #                    size factors factors,
  #                    dispersion estimates,
  #                    statistical tests ...
    DESeq()

##At this point we can have the available test results listed:
resultsNames(dds_genotype)

## A Bold Try: Model the Data With Genotype Only

##We want to see the genotype effect ...

res_genotype <- 
    results(dds_genotype, 
            name="genotype_nmrhas2_vs_creer"
            )

saveRDS(dds_genotype, file="dds_genotype.rds")

dt <- DT::datatable(
           res_genotype |> as.data.frame() |>
           dplyr::filter(!is.na(padj), padj <= 0.05) |>
           dplyr::arrange(padj),

           options=list(scrollY="500px")
     ) 

# convert numeric columns to scientific notation
dt <- dt |>
      DT::formatSignif(
           purrr::map_lgl(dt$x$data, is.numeric),
           digits=2
      )

## A Bold Try: Model the Data With Genotype Only


dt

## Why So Few Significant Genes?

##I had on purpose omitted a step which should really be done in the beginning: __`Get to know the sources of variability in your count matrix!`__ All known such sources should be either eliminated or stated in the design given to DESeq2.  Otherwise the **neglected factors will be treated as random variance** by the model and **this may mask a real effect**.  

##__`Principal Component Analysis (PCA)`__ projects the transcriptome of each sample onto a point on a two-dimensional canvas, such that the among-sample variance is maximal. 
##
##It helps to answer questions like:
##
##
##
##- Do the samples cluster as your hypotheses suggest? 
##- If not, does the cluster structure suggest __`hidden factors`__ which can be identified and included in the model? 
##- Or does __`poor clustering of expected groups`__ suggest experimental issues?

## Preparation for PCA

##The DESeq2 package provides a so-called __`"variance stabilizing transformation"`__, which reduces the increase in variance with count which is natural to count data. We will use the function to transform the input to PCA and heatmap functions.
##
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


## What Drives the Structure?

The unlabelled dataset shows some very clear structures:
plotPCA2(vsd_blind,
         intgroup = "group",
         ntop=500, 
         PCs=c(x=1,y=2))

##Genotype is more or less unrelated to this global structure!
plotPCA2(vsd_blind,
         intgroup=c("genotype"),
         color=genotype,
         ntop=500, 
         PCs=c(x=1,y=2), )

## What Drives the Structure?

plotPCA2(vsd_blind,
         intgroup="age",
         ntop=500, 
         PCs=c(x=1,y=2), 
         color = age
         )
plotPCA2(vsd_blind,
         intgroup="sex",
         ntop=500, 
         PCs=c(x=1,y=2), 
         color = sex
         )
plotPCA2(vsd_blind,
         intgroup = "tissue",
         ntop=500, 
         PCs=c(x=1,y=2),
         color = tissue
         )

## Heatmap of the 500 Most Highly Expressed Genes

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

##- With the exception of spleen and WAT, all tissues have strongly distinct expression patterns
##- Both levels of the genotype factor are compatible with any of these global patterns


## Boxplots of the Per-Sample Count Distributions 

##We want to plot the raw and normalized count distributions of each sample side by side, using `ggplot()`. For this we have to "melt" the log10 transformed matrices:


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

## Then plot it:

ggplot(df,
       aes(y=sample,x=count,fill=group)
       ) +
    geom_boxplot(show.legend=FALSE) +
    facet_wrap(vars(mode),ncol=2) +
    theme_classic() +
    theme(axis.text=element_text(size=2)) 

## Boxplots of the Per-Sample Count Distributions 


## Boxplots of the Per-Sample Count Distributions 

##- While raw counts may differ in shape and location between samples, these **differences should go away after a successful** [DESeq2 "median of ratios" normalization](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)**.
##
## ** __`Failure to normalize indicates an issue!`__. 
##
##- Here, the reason is very likely the fact that **the samples come from different tissues, with very different gene expression signatures**. DESeq2's approach to use the __`geometric expression mean over all samples of each gene as a pseudo-reference count may not be justified in this case`__. 
##

## Boxplots of the Per-Sample Count Distributions 

##To get an impression of the source of the distribution differences, we can look at the **ordered mean counts over the replicates of each level configuration**:


cnt <- counts(dds_genotype,normalized=TRUE)
medians_log10 <- log10(colMedians(cnt |> as.matrix())+1)

tapply(medians_log10, # what to group
       sub("_\\d+$","",
           names(medians_log10)
          ), # by what to group: sample name  without replicate ID
       mean  # function to apply to grouped values
       ) |> 
  sort(decreasing=TRUE) # sort the result


##- __`intestine`__ has very consistently the __`highest medians`__, followed by **WAT**
##- __`liver`__ tends to be at the __`low end`__


##- It may be advisable to split the data matrix at least into "intestine" and the remaining samples
##- __`I am not going to do this in this presentation -- YOU will do it tomorrow in the project session :)`__
##- Here I will continue with the full matrix, and see whether a better model helps to recover more of the genotype effect

## But Before We Go On ... 

## ... Will Ali Introduce You to the Math Behind DESeq2


Model formulae like  __`~ 1 + genotype`__ provide  a **high level language**

##- for the **experimenter** to __`conceptually express a hypothesis`__ about his or her data
##- for **modeling software** to know the __`structure of the tests to compute`__
##- **Independent of the actual mathematics of the tests !**

## An Additive Model of Gene Expression Including All Factors

##Including all factors in the model helps to properly capture their effects on the observed read counts. 
##This in turn may expose more of a significant, although small, effects of the genotype factor.

dds_additive <-
  # High level step 1: set up the object structure
    DESeqDataSetFromMatrix(countData = M,
                           colData = colData,
                           design= ~ age + sex + tissue + genotype 
                           ) |>
    DESeq()

resultsNames(dds_additive)
saveRDS(dds_additive, file="dds_additive.rds")


## Let's Have a Look at the Genotype Effect in this Model
res_additive <- 
    results(dds_additive, name="genotype_nmrhas2_vs_creer")

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
      DT::formatSignif(purrr::map_lgl(dt$x$data, is.numeric),
                       digits=2
      )


dt

## Yet Another Alternative: Group Contrasts

##Is there a way to compare groups of samples which differ only in the level of genotype -- in the context of the full dataset, taking the levels of all experimental factors into account?
##Breaking this sentence up into its components, we observe
##        
##- we need new factors in our colData which describe co-occurring levels of different factors in the same sample
##- we need a model in which we can compare levels of such compound factors
## Yet Another Alternative: Group Contrasts

##First we define two different compound factors:
 
##options(width=150) ## increase output line length 

  
colData <-
  colData |> 
  tidyr::unite(col=background, # name of the new column
               sep=".",
               remove=FALSE, # do NOT remove the constituent columns!
               age, sex, tissue # anything which is not a parameter is interpreted as the name of a column to unite
               ) 
head(colData,n=1)

colData <-
  colData |> 
  tidyr::unite(col=group, 
               sep=".",
               remove=FALSE,
               age, sex, tissue, genotype # genotype last for a reason!
               ) 
head(colData,n=1)

##Next we express our desired pairwise comparisons in a table:

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
nrow(contrasts)



##Observe that both the "A" and the "B" column in the table contain valid levels of the "group" factor. **For each row in the table we would like to compare samples with the "A" level combination to samples with the "B" level combination** (hence effectively comparing the effect of genotype in a fixed "A" background).

##The solution consists of two steps:


##- Run __`DESeq()`__ with a __`0 + group`__ model.  P-values from this model refer to whether or not the mean read count at a  given "group" level is significantly different from zero. Importantly the calculations take the count distribution over all levels into account. 
##- Then when passing the DESeq() result to the __`results()`__ function, we can use its __`"contrast" argument` to specifically compare a pair of "group" levels`__ listed in a table row. Although **each such comparison** involves only two levels, it **will inherit information on the full distribution from the DESeq() step**.

## Yet Another Alternative: Group Contrasts

##__`Step 1:`__

dds_group <-
  DESeqDataSetFromMatrix(countData = M,
                         colData = colData,
                         design= ~ 0 + group) |>
  DESeq()


dds_group <- readRDS("dds_group.rds")

##__`Step 2:`__

res_contrasts <- readRDS("res_contrasts.rds") 

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


##The "old.female.intestine" background has more (3692) up-regulated genotype-responsive genes than any other background. In addition it has the second highest number (1922) of down-regulated genotype-responsive genes. We will have a quick look at enriched Gene Ontology terms in these two sets of genes, to see whether they contain interesting topics. [Tomorrow we will talk in depth about Functional Enrichment.]

require("org.Mm.eg.db")

n <- "old.female.intestine"
genes <- rownames(res_contrasts[[n]])

# extract a named vector of adjusted p-values:
x <- res_contrasts[[n]] |> pull(padj) |>  
     setNames(genes)
# same for log2FoldChange:
lfc <- res_contrasts[[n]] |> pull(log2FoldChange) |>  
       setNames(genes)

L <- names(x[!is.na(x)])
l <- names(x[!is.na(x) & (x <= 0.05) & (lfc >= 1)])

clusterProfiler::enrichGO(gene = l,
                          universe = L,
                          OrgDb = org.Mm.eg.db,
                          keyType="SYMBOL",
                          ont="BP"
                  ) |> clusterProfiler::dotplot()

## Yet Another Alternative: Group Contrasts

L <- names(x[!is.na(x)])
l <- names(x[!is.na(x) & (x <= 0.05) & (lfc <= -1)])

clusterProfiler::enrichGO(gene = l,
                          universe = L,
                          OrgDb = org.Mm.eg.db,
                          keyType="SYMBOL",
                          ont="BP"
                  ) |> clusterProfiler::dotplot()

## Yet Another Alternative: Group Contrasts

n <- "old.female.liver"
genes <- rownames(res_contrasts[[n]])

# extract a named vector of adjusted p-values:
x <- res_contrasts[[n]] |> pull(padj) |>  
     setNames(genes)
# same for log2FoldChange:
lfc <- res_contrasts[[n]] |> pull(log2FoldChange) |>  
       setNames(genes)

L <- names(x[!is.na(x)])
l <- names(x[!is.na(x) & (x <= 0.05) & (lfc >= 1)])

clusterProfiler::enrichGO(gene = l,
                          universe = L,
                          OrgDb = org.Mm.eg.db,
                          keyType="SYMBOL",
                          ont="BP"
                  ) |> clusterProfiler::dotplot()

L <- names(x[!is.na(x)])
l <- names(x[!is.na(x) & (x <= 0.05) & (lfc <= -1)])

clusterProfiler::enrichGO(gene = l,
                          universe = L,
                          OrgDb = org.Mm.eg.db,
                          keyType="SYMBOL",
                          ont="BP"
                  ) |> clusterProfiler::dotplot()
