# Automated differential expression analysis for Illumina beadarrays (HT12 V4 arrays).
# Uses idat files as input.
# Author: Anand Mayakonda
# Usage: 
# source("AnalyzeBead.R")
# topTable = beadAnalyze(idats = c("file1.idat","file2.idat","file3.idat","file4.idat"),names = c("control1","control2","treated1","treated2"),
# condition = c("control","control","treated","treated"),ref.condition = "treated", fdr = 0.05, plotPCA = T)
#
# Arguments:
# idats = vector of idat files to be analyzed
# names = sample names
# condition = condition status of each sample
# ref.condition = reference condition to be used
# fdr = fdr cutoff to annotate Differentially Expressed Genes.
# plotPCA = if TRUE performs PCA and plots first two PCs
#
# returns toptable (all genes with expression and annotation, ordered according to FDR)

##########################################################################################################

#Required packages
pkgs = c('beadarray', 'illuminaHumanv4.db', 'limma', 'ggplot2', 'ggrepel')
sapply(X = pkgs, require, character.only = T)


#function to automate differential expression (from limmaDE in beadarray).
limma.autoDE = function (summaryData, SampleGroup, DesignMatrix = NULL, makeWts = TRUE,  ...) 
{
  if (is.null(SampleGroup)) 
    stop("You must define a SampleGroup for the differential expression\n")
  if (SampleGroup %in% colnames(pData(summaryData))) 
    SampleGroup <- pData(summaryData)[, SampleGroup]
  else {
    print(paste(colnames(pData(summaryData)), collapse = " "))
    stop("The SampleGroup argument must be a column name in the phenoData slot. See above for list of valid strings")
  }
  design <- model.matrix(~0 + as.factor(SampleGroup))
  colnames(design) <- as.character(levels(as.factor(SampleGroup)))
  contrast <- vector()
  for (a in 1:length(levels(as.factor(SampleGroup)))) {
    for (b in 1:length(levels(as.factor(SampleGroup)))) {
      if (a != b) {
        if (a < b) {
          contrast[length(contrast) + 1] <- paste(levels(as.factor(SampleGroup))[a], 
                                                  levels(as.factor(SampleGroup))[b], sep = "-")
        }
      }
    }
  }
  if (makeWts) {
    wts <- arrayWeights(exprs(summaryData), design)
    message("Calculating array weights")
    message("Array weights")
    wts
    fit <- lmFit(exprs(summaryData), design, weights = wts)
  }
  else fit <- lmFit(exprs(summaryData), design)
  cnt <- paste(colnames(design)[1], colnames(design)[2], sep = "-")
  cMat <- makeContrasts(contrasts = contrast, levels = design)
  fit2 <- contrasts.fit(fit, cMat)
  efit <- eBayes(fit2)
  #   deResults <- new("limmaResults", LogFC = efit$coef, PValue = efit$p.value, 
  #                    LogOdds = efit$lods, featureData = summaryData@featureData, 
  #                    phenoData = summaryData@phenoData)
  #de.df = topTable(fit2, coef=1, adjust="BH",number = 'all')
  #de.df = data.frame(LogFC = efit$coef, PValue = efit$p.value, LogOdds = efit$lods)
  return(efit)
}

#function to perform and plot PCA results
pca_plot = function(mat, pdata, label = F){
  mat.prcomp = prcomp(x = t(mat))
  pca.dat = mat.prcomp$x[,1:2]
  pca.dat = cbind(pca.dat, pdata[rownames(pca.dat),])
  if(label){
    p = ggplot(data = pca.dat, aes(x = PC1, y = PC2, color = sampleFac, label = sectionNames))+geom_point()+geom_text_repel()+theme(legend.position = 'bottom')  
  }else{
    p = ggplot(data = pca.dat, aes(x = PC1, y = PC2, color = sampleFac))+geom_point()+theme(legend.position = 'bottom')
  }
  print(p)
}

#main function to run.
beadAnalyze = function(idats, names, condition, ref.condition, fdr = 0.05, plotPCA = F){
  
  bead.batch = readIdatFiles(idatFiles = idats) #read idat files
  pData(bead.batch)[,1] = names #change default name to user provided names
  rownames(pData(bead.batch)) = names
  pData(bead.batch)[,"sampleFac"] = condition #add conditions
  pData(bead.batch)[,"sampleFac"] = as.factor(pData(bead.batch)[,"sampleFac"]) #change into factors
  pData(bead.batch)[,2] = relevel(pData(bead.batch)[,2],ref = ref.condition) #relevel the factors
  
  #normalize (available options: options are "quantile", "qspline", "vsn", "rankInvariant", "median" and "none"). Change accordingly if necessary.
  bead.eset = normaliseIllumina(BSData = bead.batch,method = "quantile",transform = "log2") 
  
  eset.exprs = as.data.frame(exprs(bead.eset)) #get exprs
  names(eset.exprs) = names
  
  is.na(eset.exprs) = do.call(cbind,lapply(eset.exprs, is.infinite)) #remove inf values if any.
  eset.exprs = as.matrix(eset.exprs[complete.cases(eset.exprs),])
  
  phenoData = new(Class = 'AnnotatedDataFrame',data = pData(bead.batch)) # create new pData
  
  eset = ExpressionSet(assayData = as.matrix(eset.exprs),phenoData = phenoData,annotation = 'Humanv4') #create new expressionSet object
  eset = addFeatureData(eset,toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE", "GENOMICLOCATION")) #add other features from IlluminaV4 pacakge.
  exprs.df = cbind(exprs(eset),as(eset@featureData,Class = 'data.frame'))
  exprs.df = exprs.df[,-grep(pattern = 'Row.names',x = colnames(exprs.df))]
  
  if(length(pData(eset)[,1]) == 2){
    cat('no replicates available; can not perform differntial expression analysis \n returning just expression table')
    return(exprs.df)
  } else{
    eset.de = limma.autoDE(summaryData = eset,SampleGroup = 'sampleFac') #run differnetial analysis.
    tt = topTable(eset.de,number = 'all')  #get toptable.
    res = merge(exprs.df,tt,by = 'row.names') #merge everything.
    res = res[order(res$adj.P.Val),] #order according to fdr.
    rownames(res) = res$Row.names
    res = res[,-1]
    
    #get up and down probes
    downProbes = rownames(res[which(res$adj.P.Val < fdr & res$logFC < 0),])
    upProbes = rownames(res[which(res$adj.P.Val < fdr & res$logFC > 0),])
    
    #annotate deg status
    res$down = as.integer(rownames(res) %in% downProbes)
    res$up = as.integer(rownames(res) %in% upProbes)
    res$up = factor(x = res$up, levels = c(1,0), labels = c('up', 'none'))
    res$down = factor(x = res$down, levels = c(1,0), labels = c('down', 'none'))
    res$deg = interaction(res$up, res$down)
    res$deg = as.character(factor(x = res$deg, levels = c("up.down"  ,"none.down","up.none"  ,"none.none"), labels = c('none', 'down', 'up', 'none')))
    res = subset(res, select = -c(up, down))
    
    #Perform Principal Component Analysis (PCA) and plot first two Principal Components
    if(plotPCA){
      pca_plot(mat = res[,names], pdata = pData(eset), label = T)  
    }
    
    #return(list(results = res, eset = eset))  
    return(res)
  }
}