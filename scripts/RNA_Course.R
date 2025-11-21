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

###Setting Up the Input for DESeq2
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


###DESeq for more factors, may expose more of a significant, although small, effects of the genotype factor
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
