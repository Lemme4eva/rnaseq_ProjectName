source("R/volcano.R")


column_boxplots <-
    function(dds) {
        require("reshape2")
        
        ## melt the input matrix
        df <-
            rbind(reshape2::melt(log10(counts(dds,normalized=FALSE)+1)) |>
                  mutate(mode="raw"),
                  
                  reshape2::melt(log10(counts(dds,normalized=TRUE )+1)) |>
                  mutate(mode="normalized")
                  ) |>
            rename(gene=Var1, sample=Var2, count=value) |>
            mutate(group=sub("_\\d+$","",sample)) |>
            mutate(mode = factor(mode,
                             levels=c("raw","normalized")
                             )
               )
        
        
        ggplot(df,
               aes(y=sample,x=count,fill=group)
               ) +
            geom_boxplot(show.legend=FALSE) +
            facet_wrap(vars(mode),ncol=2) +
            theme_classic()


    }

subset_counts <-
    function(dds,
             level_list=NULL,
             include=FALSE,
             refit=FALSE) {

        if(is.null(level_list)) {
            return(dds)
        } else if (!all(names(level_list) %in% colnames(colData(dds)))) {
            stop("Not all names of level_list are in colnames(colData(dds))")
        } else if (!all(sapply(names(level_list),
                               function(n) all(level_list[[n]] %in% colData(dds)[,n])
                               )
                        )
                   ) {
            stop("Not entries of level_list are in colData(dds)!")
        }
        keep <- Reduce(`|`,
                       lapply(names(level_list),
                              function(n)colData(dds)[,n] %in% level_list[[n]]
                              )
                       )
        if(!include) keep  <- !keep

        dds_sub <- dds[,keep]
        
        if(refit) {
            use_design <- design(dds)
            
            if(include) {
                ## if SINGLE levels are extracted (silently ASSUMED!),
                ## then the corresponding factors become meaningless
                ## in the design: remove them!

                ## -------------------------------------------------
                ## NOTE: This simple solution works only 
                ## for all-"+" or all-"*" designs!!
                ## -------------------------------------------------
                
                f <- as.character(design(dds)) ## a vector of length 2;
                                               ## [1] is "~",
                                               ## [2[ is the expression
                ## variables in level_list should be dropped
                do_drop <- intersect(strsplit(f[2],
                                              "\\s*[+*]\\s*")[[1]],
                                     names(level_list))
                
                ## first simply delete each variable:
                for(v in do_drop) f[2] <- sub(v,"",f[2])

                ## this will create an invalid formula --
                ## correct by removing adjacent duplicate components:
                x <- strsplit(f[2],"\\s+")[[1]] 
                keep <- x[c(TRUE,!(x[-1] == x[-length(x)]))]

                use_design <- as.formula(paste(f[1],
                                               paste(keep, collapse=" ")
                                               )
                                         )
                                                          
            }
            
            DESeqDataSetFromMatrix(countData = counts(dds_sub,
                                                      normalized=FALSE),
                                   colData = colData(dds_sub),
                                   design  = use_design
                                   ) |>
            DESeq()
            
        } else {
            dds_sub
        }
    }

volcano_plot <-
    function(res, check_label_overlap=TRUE,...) {

        volcano(res |> as.data.frame() |> tibble::rownames_to_column("gene"),
                use_cols=c("gene","log2FoldChange","padj"),
                lfc_thrsh=1,
                p_thrsh=0.05,
                highlight.alpha=1.0,#0.75,
                label.show=TRUE,
                check_label_overlap=check_label_overlap,
                ...)
    }


runORA <-
    function(genes, universe,OrgDb,keyType,
             
             ontologies   = c("BP","MF","CC"),
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             qvalueCutoff = 02,
             ...
             ) {


        stopifnot(all(ontologies %in% c("BP","MF","CC")))
        
        # return a list of results, one per ontology
        sapply(ontologies,
               function(ont) {
                   clusterProfiler::enrichGO(
                                        gene=genes,
                                        universe=universe,
                                        OrgDb = OrgDb,
                                        keyType=keyType,
                                        ont=ont)
               },simplify=FALSE)
                 
}


runGSEA <-
    function(geneList,OrgDb,keyType,
             
             ontologies   = c("BP","MF","CC"),
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             seed         = TRUE,
             by           = "fgsea",
             nPermSimple  = 10000,
             eps          = 1e-50,
             verbose      = FALSE,
             ...
             ) {


        stopifnot(all(ontologies %in% c("BP","MF","CC")))
        
        # return a list of results, one per ontology
        sapply(ontologies,
               function(ont) {
                   clusterProfiler::gseGO(
                                        geneList = geneList,
                                        OrgDb = OrgDb,
                                        keyType = keyType,
                                        
                                        ont=ont,
                                        minGSSize    = 10,
                                        maxGSSize    = 500,
                                        pvalueCutoff = 0.05,
                                        seed         = TRUE,
                                        by           = "fgsea",
                                        nPermSimple  = 10000,
                                        eps          = 1e-50,
                                        verbose      = FALSE)
               },simplify=FALSE)
                 
}
