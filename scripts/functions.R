require("dplyr")


plotPCA2 <- function (object, intgroup = "tissue",
                      ntop = 500, returnData = FALSE,
                      PCs=c(x=1,y=2),
                      pointSize=3,
                      hjust=-0.75,vjust=0.5,
                      just_print_select=FALSE,
                      print_labels=FALSE, 
                      use_coord_fixed=FALSE,
                      ...
                      )  {
  
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    if(just_print_select) {
        return(select)
    }
    
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    df <- as.data.frame(colData(object))

    intgroup <- intersect(intgroup,
                          colnames(df)
                          )
    if(length(intgroup)==0) {
        stop("No intersection of intgroup and colnames(colData(object))!")
    }
    intgroup.df <- df[,intgroup]

    if(length(intgroup) > 1) {
        group <- factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        group <- factor(df[,intgroup])
    }
    
    d <- data.frame(pca$x[, PCs[1]], pca$x[, PCs[2]], group = group, 
                    ##intgroup.df, name = colnames(object))
                    df, name = colnames(object))
    colnames(d)[1:2] <- paste("PC",PCs,sep="")
    if (returnData) {
        attr(d, "percentVar") <- percentVar[PCs]
        return(d)
    }

    p <- 
        ggplot(data = d,
               aes_string(x = paste("PC",PCs["x"],sep=""),
                          y = paste("PC",PCs["y"],sep="")
                          )
        ) + 
        scale_shape_manual(values=c(21:25)) +
        geom_point(mapping = aes(
                            ...
                   ),
                   size = pointSize, ## takes effect only if set at last, why?
                   ##shape = c(21,24),
                   stroke = 2
                   
        ) +
        xlab(paste0("PC",PCs["x"],": ", round(percentVar[PCs["x"]] * 100), "% variance",sep="")) + 
        ylab(paste0("PC",PCs["y"],": ", round(percentVar[PCs["y"]] * 100), "% variance",sep=""))


    if(print_labels) {
        p <- p +  geom_text(aes(label=rownames(d),hjust=hjust, vjust=vjust))
    }
    if(use_coord_fixed) { 
        p <- p + coord_fixed() 
    }

    p
}


sampleDistances <-
    function(vsd) {
        require("RColorBrewer")

        sampleDists <- dist(t(assay(vsd)))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- vsd$group2 ##paste(vsd$condition, vsd$type, sep="-")
        colnames(sampleDistMatrix) <- vsd$group2 ##NULL
        colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
        pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)
    }

get_orgDb <- 
  function(organism_name) {
    orgDb <- 
      function() {
        AH <- AnnotationHub::AnnotationHub(cache = Sys.getenv("ANNOTATION_HUB_CACHE"))
        OrgDB <- AnnotationHub::query(
          AH, pattern = c(organism_name, "OrgDB")
        )[[1]]
      }
    return(orgDb)
  }

filter_result <- 
  function(res,
           padj_thrsh = 0.05, ## Call a gene significant if padj < padj_thrsh
           lfc_thrsh  = 1.5,  ## Call a gene relevant 
           ##  if abs(log2FoldChange) >= lfc_thrsh
           direction = "up" 
  ) {
    
    if(direction == "up") {
      cmp <- `>=`
    } else if (direction == "down") {
      cmp <- `<=`
    } else {
      stop("unknown direction!")
    }
    
    res |> 
      filter(padj < padj_thrsh, 
             abs(log2FoldChange) >= lfc_thrsh, 
             cmp(log2FoldChange,0)
      )
  }

run_enrichGO <-
  function(genes, universe, orgdb, ontology,keytype,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.01,
           qvalueCutoff  = 0.05,
           readable      = TRUE
           )

    enrichGO(gene          = genes,
             universe      = universe,
             OrgDb         = orgdb,
             ont           = ontology,
             keyType       = keytype,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05,
             readable      = TRUE)


