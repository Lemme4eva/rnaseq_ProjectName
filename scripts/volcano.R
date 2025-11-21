# from: /home/ugoebel/CECAD/Trifunovic/Project_Kaul/Analysis/MultiOmics/MixOmics/Testing_mixOmics/Testing_R-4.3.0/R/functionsUG.R
## A derivative of MetaboAnalystR::PlotVolcano
## Original function name: my_MA_volcano 

volcano <- function(df, ## a named data.frame with >= 3 columns,
                              use_cols = c(2,4,5),
                              ## applies to result from get_tx_l2fc)(this_pair=this_pair),
                              ## using p.value (p.adj would be column 6)
                              
                              lfc_thrsh, p_thrsh,
                              
                              highlight.what = NULL, ## a list of lists,
                                                     ## with sublist names = color names,
						     ## sublist content = points to be highlighted
                              highlight.pointsize = 4, #2,##4,#6,
                              highlight.textsize = 6,
                              highlight.alpha = 0.5,
                              
                              background.color = "grey",
                              background.pointsize = 3, #2, ##3,
                              background.textsize = 3,
                              background.alpha = 0.3,

                              nonsig.pointsize = 0.1,
                              
                              label.show = FALSE, 
                              label.box  = FALSE,
                              
                              imgName = "test_",

                              plotTheme = 0,
                              format = "png", 
                              dpi = 300, width = NA,
                              y_mode="simple",
                              y_breaks=NULL,
                              y_break_scales=NULL,
                              y_break_which =NULL,
                              y_break_expand=FALSE,
                              y_break_space=0.1,
                              
                              expand_ceiling=0,
                              check_label_overlap = TRUE,
                              ...
                    ) {


    good_rows <- which(!is.na(df[,use_cols[1]]) & ## [1]= gene(!) name can be NA with protein data
                                                  ## (not if [1] is the protein ID)
                       !is.na(df[,use_cols[3]]))

    df <- data.frame(df)
    
    df <- df[good_rows, use_cols]
    
    if(!is.null(y_breaks)) {
        ## add a dummy very small dummy p-value to df to
        ## expand the range of -log10(p) upward
        ## for better plotting of labels:
        dummy_df <- data.frame("NOSHOW",
                               0,
                               10^(- (ceiling(max(
                                         -log10(df[,3])
                                                 )
                                             ) + expand_ceiling
                                     )
                                  )
                               )
        colnames(dummy_df) <- colnames(df)
        data <- rbind(dummy_df,
                      df)
    } else {
        data <- df
    }
    

    
    NOSHOW <- data[,1] == "NOSHOW"
    UP   <-   !NOSHOW & (data[,2] >=  lfc_thrsh)
    DOWN <-   !NOSHOW & (data[,2] <= -lfc_thrsh)
    P    <-   !NOSHOW & (-log10(data[,3]) >= -log10(p_thrsh))
    data$significant <- P&(UP|DOWN)
    
    if(is.null(highlight.what)) {
        ## if not otherwise stated,
        ## highlight significant points
        highlight.what <-
            list(dummy_point =
                     list(color=NA,
                          symbols=data[which(NOSHOW),1]
                          ),
                 up_regulated =
                     list(color="firebrick",
                          symbols=data[which(P & UP),1]
                          ),
                 down_regulated =
                     list(color="cornflowerblue",
                          symbols=data[which(P & DOWN),1]
                          )
                 )
    }

    data$group <- "other"
    data$color <- background.color
    data$pointsize <- background.pointsize
    data$textsize <- background.textsize
    data$highlighted  <- FALSE
    data$alpha <- background.alpha
    data$label <- NA

    ## assign group membership to genes in the order of
    ## decreasing group size
    ## -- this gives genes in small groups a chance to be seen
    size_order <-
        order(sapply(highlight.what,
                     function(x)length(x$symbols)
                     ),
              decreasing=TRUE)
    
    for(group in names(highlight.what)[size_order]) {
        i <- which(data[,1] %in% highlight.what[[group]]$symbols)
        
        if(length(i) > 0) {
            data$group[i] <- group
            data$color[i] <- highlight.what[[group]]$color
            data$pointsize[i] <- highlight.pointsize
            data$alpha[i] <- highlight.alpha
            data$highlighted[i] <- TRUE
            
            if(label.show) {
                data$label[i] <- data[i,1] ##data[1,i] ## UG Nov 17, 2025: ?? should be data[i,1] ????
                data$textsize[i] <- highlight.textsize
            }
        }
    }
    ## finally, re-set the pointsize in the non-significant area,
    ## for all groups and for "other"
    data$pointsize[which(!data$significant)] <- nonsig.pointsize
    
    ## re-order to have highlighted points last:
    data <- rbind(data %>% filter( group=="other"),
                  data %>% filter(!group=="other"))

    data$group <- as.factor(data$group)
    group_colors <- sapply(levels(data$group),
                           function(x)unique(
                                          data$color[data$group==x]
                                      )
                           ) ## ugly, but for now ...
                             ## mind possible issue with == and factor ??

    data$plotregion <- factor(ifelse(!data$significant,
                                     "non-sig",
                              ifelse(data$highlighted,
                                     "of interest", ##"gene sets",
                                     "background genes"
                                     )
                              ))

                              
    ##data$highlighted <- as.factor(data$highlighted)

    alphas <- c("background genes"=background.alpha,
                "of interest" =highlight.alpha,
                "non-sig"   =background.alpha)
    pointsizes <- c("background genes"=background.pointsize,
                    "of interest" =highlight.pointsize,
                    "non-sig"   =nonsig.pointsize)


##    imgName = paste(imgName, "_dpi", dpi, ".", format, sep = "")
##    if (is.na(width)) {
##        w <- 10
##    }
##    else if (width == 0) {
##        w <- 8
##    }
##    else {
##        w <- width
##    }
##    h <- w * 6/10
##    
##    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w, 
##       height = h, type = format, bg = "white")

    require(ggplot2)
    p <- ggplot(data = data,
                aes(x = data[, 2],
                    y = -log10(data[, 3]),
                    col = group,
                    alpha = plotregion,
                    size  = plotregion, 
                    label = label)
                ) +
        scale_color_manual(values = group_colors) +
        scale_alpha_manual(values = alphas) + 
        scale_size_manual(values = pointsizes) + 
        
        geom_vline(xintercept = c(-lfc_thrsh, lfc_thrsh),
                   linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(p_thrsh), 
                   linetype = "dashed", color = "black") +
        geom_point() +
        #geom_label() + #check_overlap = TRUE) +
        geom_text(check_overlap = check_label_overlap, color="black") +
                
                
        labs(x = "log2(FC)", y = "-log10(p)")



    if (plotTheme == 0) {
        p <- p + theme_bw() ## an alias of ggplot2::theme_grey
    }
    else if (plotTheme == 1) {
        p <- p + theme_grey()
    }
    else if (plotTheme == 2) {
        p <- p + theme_minimal()
    }
    else {
        p <- p + theme_classic()
    }
    
    if(!is.null(y_breaks)) {
        p <- p + ggbreak::scale_y_cut(
                              space=y_break_space,
                              breaks=y_breaks,
                              scales=y_break_scales,
                              which=y_break_which,
                              expand= y_break_expand ## this adds the space!
                          )
    } 

    p
##    print(p)
##    dev.off()
    
}
