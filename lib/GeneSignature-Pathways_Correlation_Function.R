

### Read in Functional Gene-Set GMT

readGMT <- function(gmtfile, normalized=F){
  gmt <- readLines(gmtfile)
  geneset <- strsplit(gmt, split="\t")
  
  Geneset_names <- vapply(geneset, function(y) y[1], character(1))
  
  # extract the all names as new name
  for (i in 1:length(Geneset_names)){
    name <- Geneset_names[i]
    name <- unlist(strsplit(name, split="_"))
    name <- name[c(1:length(name))]
    name <- paste(name, collapse = "_")
    Geneset_names[i] <- name
  }
  names(geneset) <- Geneset_names
  
  if (normalized == FALSE){
    geneset <- lapply(geneset, "[", -c(1:2))
  } else if (normalized == TRUE){
    geneset <- lapply(geneset, "[", -c(1))
  }
  
  return(geneset)
}


### normalize Gene-Set name

normalize_GeneSet_Name <- function(geneSet, geneNameAnno){
  # set progress bar
  message("提示：正在将基因统一规范为Symbol：")
  progressBar <- txtProgressBar(style=3)
  
  # 遍历每一个基因集
  for (k in 1:length(geneSet)){
    # progress bar
    setTxtProgressBar(progressBar, k/length(geneSet))
    
    one_geneSet <- geneSet[[k]]
    
    # 遍历基因集中每个基因
    for (i in 1:length(one_geneSet)){
      # 检查是否为Gene Symbol
      if (!(one_geneSet[i] %in% names(geneNameAnno))){
        # 检查是否为别称
        for (j in 1:length(geneNameAnno)){
          if (one_geneSet[i] %in% unlist(geneNameAnno[j])){
            # 更正名称为Symbol
            one_geneSet[i] <- names(geneNameAnno)[j]
          }
        }
      }
    }
    
    geneSet[[k]] <- one_geneSet
  }
  # progress bar close
  close(progressBar)
  
  message(paste("提示：完成！共转换了",as.character(k),"个基因集的基因名",sep=""))
  
  return(geneSet)
}

###############################################################################
### Score calculation functions ----
###############################################################################



get_GeneSignature_Score <- function(data, genes, pickMethod="mean", log="log2(n)", exp=FALSE, calOrder){
  
  # missNum <- length(genes) - length(rownames(data[rownames(data) %in% genes,]))
  # if (missNum != 0){
  #   missGene <- setdiff(genes, rownames(data[rownames(data) %in% genes,]))
  #   missGene <- paste(missGene, sep=",")
  #   message(paste("缺失了", missNum, "个基因：", missGene, "\n", sep=""))
  # }
  
  SignatureGenes_expr <- data[rownames(data) %in% genes,]
  
  if (log == "log2(n+1)"){
    SignatureGenes_expr <- log2(SignatureGenes_expr+1)
  } else if(log == "log2(n)"){
    SignatureGenes_expr <- log2(SignatureGenes_expr)
  }
  
  # extract expression for score
  #===========================================================================
  if (calOrder == "calFirst"){
    # calFirst
    if (pickMethod == "median"){
      Score <- apply(SignatureGenes_expr,2,function(tmp){
        median(tmp)
              #** 1: 按行操作， 2：按列操作
      })
    } else if (pickMethod == "mean"){
      Score <- apply(SignatureGenes_expr,2,function(tmp){
        mean(tmp)
      })
    }
    
  } else if (calOrder == "standFirst"){
    # scaleFirst
    Score <- apply(SignatureGenes_expr, 1, function(tmp){
      scale(tmp)
    })
    
    if (pickMethod == "median"){
      Score <- apply(Score, 1, function(tmp){
        median(tmp)
      })
      
    } else if (pickMethod == "mean"){
      Score <- apply(Score, 1, function(tmp){
        mean(tmp)
      })
    }
  }
  
  # fix methods
  #===========================================
  
  if (exp == TRUE){
    Score <- exp(Score)
  }
  
  # 原理：每个样本中取Signature Genes表达中值，并取幂
  # Calculate GeneSignature raw scores by exponentiation of the median value in the Genes
  
  # origin
  # Score <- apply(data[genes,],2,function(tmp){
  #   exp(median(tmp))
  # })

  return(Score)
}



get_Geneset_Score <-function(geneset, expr, pickMethod="mean", log=TRUE, exp=FALSE, calOrder){
  
  # precreate a dataframe for containing scores
  Geneset_Scores_DF <- data.frame(names=colnames(expr))
  Geneset_Scores_DF <- data.frame(t(Geneset_Scores_DF))
  colnames(Geneset_Scores_DF) <- Geneset_Scores_DF[1,]
  
  # i for progress bar of calculating GO-sets scores
  i <- 2
  
  message("提示：正在计算基因集的GeneSignature打分：")
  progressBar <- txtProgressBar(style=3)
  
  # loop every Term from GO&KEGG
  ###############################################################
  for (name in names(geneset)){
    # progress bar
    setTxtProgressBar(progressBar, i/length(geneset))
    
    # extract one set's genes
    oneGeneSet_Genes <- geneset[[name]]
    
    # calculate GeneSignature scores
    oneGeneSet_Scores_raw <- get_GeneSignature_Score(expr, oneGeneSet_Genes,
                                                     pickMethod = pickMethod,
                                                     log = log,
                                                     exp = exp,
                                                     calOrder = calOrder)
    # # scale scores
    if (calOrder == "calFirst"){
      oneGeneSet_Scores <- as.numeric(scale(oneGeneSet_Scores_raw))
      # get sample name
      names(oneGeneSet_Scores) <- names(oneGeneSet_Scores_raw)
    } else if (calOrder == "standFirst"){
      oneGeneSet_Scores <- oneGeneSet_Scores_raw
    }
    
    # merge
    Geneset_Scores_DF <- rbind(Geneset_Scores_DF, oneGeneSet_Scores)
    rownames(Geneset_Scores_DF)[i] <- name
    
    i <- i+1
  }
  close(progressBar)
  
  Geneset_Scores_DF <- Geneset_Scores_DF[-1,]
  setNames <- rownames(Geneset_Scores_DF)
  
  Geneset_Scores_DF <- as.data.frame(lapply(Geneset_Scores_DF, as.numeric))
  rownames(Geneset_Scores_DF) <- setNames
  
  message(paste("提示：完成！共计算了",length(rownames(Geneset_Scores_DF)),"个基因集的GeneSignature打分",sep=""))
  return(Geneset_Scores_DF)
}



getHighLow <- function(x, str = NULL) {
  
  # Stratify into High/Low groups
  
  classes <- rep("Lo", length(x))
  classes[which(x>median(x))] <- "Hi"
  
  if (!is.null(str)) classes <- paste(str, classes, sep="_")
  
  names(classes) <- names(x)
  classes <- as.factor(classes)
  lowclass <- levels(classes)[grep("Lo", levels(classes))]
  classes <- relevel(classes, ref=lowclass)
  
  return(classes)
  
}


### Differential expression analysis functions ----


runLimma <- function(data, HiLo, alpha=0.05) {
  
  reflev <- unique(HiLo)
  reflev <- as.character(reflev[grep("Lo", reflev)])
  group <- relevel(HiLo,ref=reflev)
  design <- model.matrix(~group)
  
  fit <- lmFit(data, design)
  fit <- eBayes(fit)
  DEGs <- topTable(fit, coef=colnames(design)[2], adjust="BH", number=nrow(data))
  
  DEGs$DE <- NA
  DEGs$DE[intersect(which(DEGs$logFC>0), which(DEGs$adj.P.Val<alpha))] <- "up"
  DEGs$DE[intersect(which(DEGs$logFC<0), which(DEGs$adj.P.Val<alpha))] <- "down"
  
  return(DEGs)
  
}


###############################################################################
### Plot functions ----
###############################################################################


plotViolin <- function(scoreData, filename){
  
  minp <- 10^-10
  
  pdf(filename, onefile=TRUE, width = 8, height = 8)
  
  pvals <- matrix(nrow = 3, ncol = ncol(scoreData)-1)
  rownames(pvals) <- paste0("strat", 2:4)
  colnames(pvals) <- colnames(scoreData)[-1]
  for (i in 2:length(colnames(scoreData))){
    
    for (q in 2:4){
      
      if (q <= 2){
        res <- kruskal.test(scoreData[,i] ~ getStrata(scoreData$ATG, q))
        ptext <- paste0("pval=", format.pval(res$p.value, digits = 4, eps = minp))
        plotcols <- ATGcolors
      } else {
        res <- kruskal.test(scoreData[,i] ~ getStrata(scoreData$ATG, q))
        ptext <- paste0("pval=", format.pval(res$p.value, digits = 4, eps = minp))
        plotcols <- rev(colorRampPalette(ATGcolors)(q))
        names(plotcols) <- levels(getStrata(scoreData$ATG, q))
      }
      
      plotData <- cbind(ATG = getStrata(scoreData$ATG, q), scoreData[,-1,drop=FALSE])
      print(ViolinPlot(plotData, c("ATG", colnames(scoreData)[i]), cols = plotcols, ptext))
      
      pvals[q-1,i-1] <- res$p.value
    }
    
  }
  
  dev.off()
  
  return(pvals)
}


HeatMap <- function(heatdata, titlestr, annDF, annColors, file, dendrogram = NULL){
  
  mycols <- colorRamp2(breaks = c(-2, 0, 2), colors=c("#313695", "white", "#d73027"))
  
  ha <- HeatmapAnnotation(df=annDF, col=annColors)
  
  if (!is.null(dendrogram)){
    ccol <- dendrogram
  } else {
    ccol <- TRUE
  }
  
  hm <- Heatmap(data.matrix(heatdata), 
                name="log2FC", 
                col=mycols,
                show_row_names=FALSE,
                show_column_names=FALSE,
                column_title=titlestr, 
                row_title="Differentially Expressed Genes",
                top_annotation=ha,
                cluster_rows = FALSE,
                cluster_columns = ccol)
  
  pdf(file, onefile=TRUE, width=12, height=8)
  draw(hm)
  dev.off()
  
  return(hm)
  
}