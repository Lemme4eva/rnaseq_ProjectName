plotPCA2 <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE,
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
        scale_shape_manual(values=c(21,24)) +
        geom_point(mapping = aes(
                       ...
                   ),
                   size = pointSize, ## takes effect only if set at last, why?
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
