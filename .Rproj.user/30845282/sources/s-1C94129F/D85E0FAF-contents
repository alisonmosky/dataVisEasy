####Statistics and Functions for Correlations and Expressions####

#' @importFrom graphics hist
#' @importFrom stats TukeyHSD aov as.dist cor cor.test cutree density dist hclust mad median t.test
#' @importFrom reshape2 melt
#' @importFrom dplyr summarise_all summarise_at group_by bind_rows

#' @export
AOV1way <- function(
  data.to.aov,
  category,
  pthreshold = 0.05,
  additional.report = "NONE"  ##options are "NONE", "TUKEY","AOV", or "ALL"
){

  if (("matrix" %in% class(data.to.aov)) != TRUE ) {
    data.to.aov <- as.matrix(data.to.aov)
    warning('input data converted to matrix')
  }

  ####if data is not all samples, subset annotations appropriately
  if (sum(colnames(data.to.aov) %notin% rownames(params$annotations)) != 0 ) {
    stop('colnames of input data do not match rownames of annotations, cannot link annotations to data and assign groupings for ANOVA')}

  temp.annotations <- params$annotations[match(colnames(data.to.aov), rownames(params$annotations)),]

  groupings <- droplevels(as.factor(temp.annotations[,category]))


  ##remove rows where there are no observations for one category


  dat.check <- t(data.to.aov); dat.check <- reshape2::melt(dat.check); dat.check$groupings <- groupings
  present <- dplyr::summarise_all(dplyr::group_by(dat.check,Var2, groupings), list(present=function(x)(sum(!is.na(x)))))
  incomplete <- unique(as.character(present$Var2[which(present$value_present == 0)]))


  present.incomplete <- present[which(present$Var2 %in% incomplete),] %>% group_by(Var2)
  levels.present <- present.incomplete %>% dplyr::summarise_at("value_present",function(x)(sum(x != 0)))
  take.out <- unique(as.character(levels.present$Var2[which(levels.present$value_present <= 1)]))

  if (length(take.out != 0)) {
    data.to.aov <- data.to.aov[-which(rownames(data.to.aov) %in% take.out),]
    warning(paste(paste(take.out, collapse = ", "),"do not have two or more groups with no non-missing arguments and have been removed from the AOV analysis "))

  }

  incomplete.keep <- incomplete[incomplete %notin% take.out]
  if (length(incomplete.keep) != 0 ) {
    warning(paste(paste(incomplete.keep, collapse = ", "), "have at least one group with no non-missing arguments"))
  }



  aov.all <- apply(data.to.aov, 1, function(x)(summary(aov(x~groupings))))

  aov.results <- data.frame(FVal=unlist(lapply(aov.all,function(x)((x[[1]]$`F value`[1])))),
                            pVal=unlist(lapply(aov.all,function(x)((x[[1]]$`Pr(>F)`[1])))),
                            row.names = names(aov.all))

  sig.genes <- rownames(aov.results)[which(aov.results$pVal <= pthreshold)]
  nonsig.genes <- rownames(aov.results)[which(aov.results$pVal > pthreshold)]

  sig.set <- data.to.aov[which(rownames(data.to.aov) %in% sig.genes),]

  tukey.all <-  apply(sig.set,1,function(x)(TukeyHSD(aov(x~groupings))))

  tukey.pvals <- data.frame(Reduce(dplyr::bind_rows, lapply(tukey.all,function(x)(x$groupings[,4]))), row.names = names(tukey.all)); colnames(tukey.pvals) <- gsub("\\.","-", colnames(tukey.pvals))
  tukey.diffs <- data.frame(Reduce(dplyr::bind_rows, lapply(tukey.all,function(x)(x$groupings[,1]))), row.names = names(tukey.all)); colnames(tukey.diffs) <- gsub("\\.","-", colnames(tukey.diffs))


  if (toupper(additional.report) == "ALL") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                'Sig.Genes' = sig.genes,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.output' = tukey.all,
                'Tukey.pVals' = tukey.pvals),
           'Tukey.diffs' = tukey.diffs)
  }

  if (toupper(additional.report) == "TUKEY") {
    return(list('AOV.Results' = aov.results,
                'Sig.Genes' = sig.genes,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.output' = tukey.all,
                'Tukey.pVals' = tukey.pvals,
                'Tukey.diffs' = tukey.diffs))
  }

  if (toupper(additional.report) == "AOV") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                'Sig.Genes' = sig.genes,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.pVals' = tukey.pvals,
                'Tukey.diffs' = tukey.diffs))
  }

  if (toupper(additional.report) == "NONE") {
    return(list('AOV.Results' = aov.results,
                'Sig.Genes' = sig.genes,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.pVals' = tukey.pvals,
                'Tukey.diffs' = tukey.diffs))
  }

}



#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom plyr join_all
#' @export
AOV2way <- function(
  data.to.aov,
  category1,
  category2,
  pthreshold = 0.05,
  additional.report = "NONE"  ##options are "NONE", "TUKEY","AOV", or "ALL"
){

  if (("matrix" %in% class(data.to.aov)) != TRUE ) {
    data.to.aov <- as.matrix(data.to.aov)
    warning('input data converted to matrix')
  }

  ####if data is not all samples, subset annotations appropriately
  if (sum(colnames(data.to.aov) %notin% rownames(params$annotations)) != 0 ) {
    stop('colnames of input data do not match rownames of annotations, cannot link annotations to data and assign groupings for ANOVA')}

  temp.annotations <- params$annotations[match(colnames(data.to.aov), rownames(params$annotations)),]

  groupings1 <- droplevels(as.factor(temp.annotations[,category1]))
  groupings2 <- droplevels(as.factor(temp.annotations[,category2]))


  ##remove rows where there are no observations for either category


  dat.check <- t(data.to.aov); dat.check <- melt(dat.check); dat.check$grouping1 <- groupings1; dat.check$grouping2 <- groupings2
  present1 <- dplyr::summarise_all(dplyr::group_by(dat.check,Var2, grouping1), list(present=function(x)(sum(!is.na(x)))))
  present2 <- dplyr::summarise_all(dplyr::group_by(dat.check,Var2, grouping2), list(present=function(x)(sum(!is.na(x)))))

  incomplete1 <- unique(as.character(present1$Var2[which(present1$value_present == 0)]))
  present.incomplete1 <- present1[which(present1$Var2 %in% incomplete1),] %>% group_by(Var2)
  levels.present1 <- present.incomplete1 %>% dplyr::summarise_at("value_present",function(x)(sum(x != 0)))
  take.out1 <- unique(as.character(levels.present1$Var2[which(levels.present1$value_present <= 1)]))

  incomplete2 <- unique(as.character(present2$Var2[which(present2$value_present == 0)]))
  present.incomplete2 <- present2[which(present2$Var2 %in% incomplete2),] %>% group_by(Var2)
  levels.present2 <- present.incomplete2 %>% dplyr::summarise_at("value_present",function(x)(sum(x != 0)))
  take.out2 <- unique(as.character(levels.present2$Var2[which(levels.present2$value_present <= 1)]))

  if (length(unique(c(take.out1,take.out2)) != 0)) {
    data.to.aov <- data.to.aov[-which(rownames(data.to.aov) %in% unique(c(take.out1,take.out2))),]
    warning(paste(paste(unique(c(take.out1,take.out2)), collapse = ", "),"do not have two or more groups with no non-missing arguments and have been removed from the AOV analysis "))

  }

  incomplete.keep <- unique(c(incomplete1,incomplete2))[unique(c(incomplete1,incomplete2)) %notin% unique(c(take.out1,take.out2))]
  if (length(incomplete.keep) != 0 ) {
    warning(paste(paste(incomplete.keep, collapse = ", "), "have at least one group with no non-missing arguments"))
  }


  ###check if all effects are estimable
  aov.check <- aov(data.to.aov[1,]~groupings1 + groupings2 + groupings1:groupings2)
  if (aov.check$rank < length(aov.check$coefficients)) {
    stop('Some effects not estimable. This is likely due to overlap in the two categories provided. Consider using AOV1way with one of the categories provided instead')
  }

  aov.all <- apply(data.to.aov, 1, function(x)(summary(aov(x~groupings1 + groupings2 + groupings1:groupings2))))


  aov.results <- t(data.frame((lapply(aov.all, function(x)(unlist(x[[1]]))))))[,c(13:15,17:19)];


  colnames(aov.results) <- c(paste0("FVal-",category1),paste0("FVal-",category2),paste0("FVal-",category1,":",category2),
                             paste0("pVal-",category1),paste0("pVal-",category2),paste0("pVal-",category1,":",category2))



  category1.sig <- rownames(aov.results)[which(aov.results[,4] <= pthreshold)]
  category2.sig <- rownames(aov.results)[which(aov.results[,5] <= pthreshold)]
  interaction.sig <- rownames(aov.results)[which(aov.results[,6] <= pthreshold)]
  any.sig <- unique(c(category1.sig, category2.sig, interaction.sig))
  nonsig.genes <- rownames(data.to.aov)[rownames(data.to.aov) %notin% any.sig]

  sig.set <- data.to.aov[which(rownames(data.to.aov) %in% any.sig),]

  tukey.all <- apply(data.to.aov, 1, function(x)(TukeyHSD(aov(x~groupings1 + groupings2 + groupings1:groupings2))))

  tukey.pvals1 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% category1.sig)],function(x)(x$groupings1[,4, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% category1.sig)]); colnames(tukey.pvals1) <- rownames(tukey.all[[1]]$groupings1) # gsub("\\.","-", colnames(tukey.pvals1))
  tukey.pvals2 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% category2.sig)],function(x)(x$groupings2[,4, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% category2.sig)]); colnames(tukey.pvals2) <- rownames(tukey.all[[1]]$groupings2) # gsub("\\.","-", colnames(tukey.pvals1))
  tukey.pvals3 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% interaction.sig)],function(x)(x$`groupings1:groupings2`[,4, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% interaction.sig)]); colnames(tukey.pvals3) <- rownames(tukey.all[[1]]$`groupings1:groupings2`) # gsub("\\.","-", colnames(tukey.pvals1))

  tukey.diffs1 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% category1.sig)],function(x)(x$groupings1[,1, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% category1.sig)]); colnames(tukey.diffs1) <- rownames(tukey.all[[1]]$groupings1) # gsub("\\.","-", colnames(tukey.diffs1))
  tukey.diffs2 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% category2.sig)],function(x)(x$groupings2[,1, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% category2.sig)]); colnames(tukey.diffs2) <- rownames(tukey.all[[1]]$groupings2) # gsub("\\.","-", colnames(tukey.diffs1))
  tukey.diffs3 <- data.frame(lapply(lapply(lapply(tukey.all[which(names(tukey.all) %in% interaction.sig)],function(x)(x$`groupings1:groupings2`[,1, drop = FALSE])), as.data.frame ), tibble::rownames_to_column) %>% plyr::join_all(by = "rowname") %>% tibble::column_to_rownames("rowname") %>%  t(), row.names = names(tukey.all)[which(names(tukey.all) %in% interaction.sig)]); colnames(tukey.diffs3) <- rownames(tukey.all[[1]]$`groupings1:groupings2`) # gsub("\\.","-", colnames(tukey.diffs1))


  if (toupper(additional.report) == "ALL") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                "Category1_Sig.Genes" = category1.sig,
                "Category2_Sig.Genes" = category2.sig,
                "Interaction_Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.output' = tukey.all,
                'Category1_Tukey.pVals' = tukey.pvals1,
                'Category2_Tukey.pVals' = tukey.pvals2,
                'Interaction_Tukey.pVals' = tukey.pvals3,
                'Category1_Tukey.diffs' = tukey.diffs1,
                'Category2_Tukey.diffs' = tukey.diffs2,
                'Interaction_Tukey.diffs' = tukey.diffs3
    ))
  }

  if (toupper(additional.report) == "AOV") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                "Category1_Sig.Genes" = category1.sig,
                "Category2_Sig.Genes" = category2.sig,
                "Interaction_Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1_Tukey.pVals' = tukey.pvals1,
                'Category2_Tukey.pVals' = tukey.pvals2,
                'Interaction_Tukey.pVals' = tukey.pvals3,
                'Category1_Tukey.diffs' = tukey.diffs1,
                'Category2_Tukey.diffs' = tukey.diffs2,
                'Interaction_Tukey.diffs' = tukey.diffs3
    ))
  }

  if (toupper(additional.report) == "TUKEY") {
    return(list('AOV.Results' = aov.results,
                "Category1_Sig.Genes" = category1.sig,
                "Category2_Sig.Genes" = category2.sig,
                "Interaction_Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1_Tukey.pVals' = tukey.pvals1,
                'Category2_Tukey.pVals' = tukey.pvals2,
                'Interaction_Tukey.pVals' = tukey.pvals3,
                'Category1_Tukey.diffs' = tukey.diffs1,
                'Category2_Tukey.diffs' = tukey.diffs2,
                'Interaction_Tukey.diffs' = tukey.diffs3
    ))
  }

  if (toupper(additional.report) == "NONE") {
    return(list('AOV.Results' = aov.results,
                "Category1_Sig.Genes" = category1.sig,
                "Category2_Sig.Genes" = category2.sig,
                "Interaction_Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1_Tukey.pVals' = tukey.pvals1,
                'Category2_Tukey.pVals' = tukey.pvals2,
                'Interaction_Tukey.pVals' = tukey.pvals3,
                'Category1_Tukey.diffs' = tukey.diffs1,
                'Category2_Tukey.diffs' = tukey.diffs2,
                'Interaction_Tukey.diffs' = tukey.diffs3
    ))
  }
}


#' @importFrom pcaMethods pca scores loadings
#' @export
myPCA <- function(
  data,
  to.pca = "samples",
  nPcs= 3,
  color.by = "blue", ##vector same length or in annotations, annot_cols
  custom.color.vec = FALSE,
  PCs.to.plot = c("PC1","PC2"),
  legend.position = "right",
  main = NULL,
  point.size =5,
  transparency = 1,
  percent.mad =0.5,
  return.ggplot.input = FALSE,
  return.loadings = FALSE
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  if (to.pca == "samples") {

    temp.annotations <- params$annotations


    pca <- pcaMethods::pca(t(data), nPcs = nPcs)
    pca.scrs <- pcaMethods::scores(pca)
    pca.ldgs <- pcaMethods::loadings(pca)



    if (color.by %in% rownames(data) | sum(custom.color.vec != FALSE) > 0) {
      pca.data <- data.frame(pca.scrs, Samples = colnames(data))
      if (color.by %in% rownames(data)) {
        genedat<- data[which(rownames(data)==color.by),]
        cols <- myColorRamp5(params$expression_gradient.colors,genedat, percent.mad = percent.mad)
      } else{ cols <- custom.color.vec}

      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols, Samples = Samples))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
        scale_fill_identity() +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)

    }else{


      if (color.by %in% colnames(temp.annotations)) {
        if (sum(!is.na(temp.annotations)) != 0) {
          if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
            stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')}

          temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]
          #temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),]
          pca.data <- data.frame(pca.scrs,temp.annotations,Samples = colnames(data))
        }

        if (color.by %in% names(params$annot_cols)) {
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
        }else{
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- scales::hue_pal()(length(levels(cols)))
        }

      } else{ cols <- color.by; colors <- color.by; pca.data <- data.frame(pca.scrs, Samples = colnames(data)); legend.position = "none"}


      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols, Samples = Samples))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
        scale_fill_manual(values=colors) +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2]), fill=color.by) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
    }
  }

  if (to.pca == "genes") {

    temp.annotations.genes <- params$annotations.genes

    pca <- pcaMethods::pca((data), nPcs = nPcs)
    pca.scrs <- pcaMethods::scores(pca)
    pca.ldgs <- pcaMethods::loadings(pca)



    if (sum(custom.color.vec != FALSE) > 0) {
      pca.data <- data.frame(pca.scrs, Genes = rownames(data))
      cols <- custom.color.vec

      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols, Genes = Genes))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
        scale_fill_identity() +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)

    }else{


      if (color.by %in% colnames(temp.annotations.genes)) {
        if (sum(!is.na(temp.annotations.genes)) != 0) {
          if (sum(rownames(data) %notin% rownames(temp.annotations.genes)) != 0 ) {
            stop('rownames of input data do not match rownames of annotations, cannot link annotations to data')}

          temp.annotations.genes <- temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),, drop = FALSE]
          pca.data <- data.frame(pca.scrs,temp.annotations.genes, Genes = rownames(data))
          #temp.annotations.genes <-temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),]
        }
        if (color.by %in% names(params$annot_cols)) {
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
        }else{
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- scales::hue_pal()(length(levels(cols)))
        }

      } else{ cols <- color.by; colors <- color.by; pca.data <- data.frame(pca.scrs, Genes = rownames(data)); legend.position = "none"}


      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols, Genes = Genes))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
        scale_fill_manual(values=colors) +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2]), fill=color.by) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
    }
  }
  if (return.loadings == TRUE) {return(pca.ldgs)}
  if (return.ggplot.input == TRUE) {return(pca.data)}
  return(p)
}




#' @export
PTM <- function(  #pavlidis template matching, correlation to a chosen template
  data, ##dataset, genes in the rows
  match.template, ##from params$annotations or gene, can also be gene or sample
  annotation.level.set.high, ##level(s) of annotations vector to set high
  custom.template = NA,  ##must be same length as data
  Find.Match.For = "genes", ##default to compare to a gene or sample metadata template, can change to "samples"
  cutoff = 0.05, ##default of 0.05, can be changed, assumes pval
  cut.by = "pvals", ##pvals will return lower than cutoff with a positive correlation, can also use rvals, returns values greater than cutoff
  method = "pearson", ##default of pearson
  NA.handling = "pairwise.complete.obs", ##default of pairwise complete obs, can be changed
  return.vals = FALSE
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  if (tolower(Find.Match.For) == "genes") {   ##length of samples

    temp.annotations <- params$annotations



    if(is.na(custom.template)==TRUE){

      if (((match.template %in% colnames(temp.annotations)) == FALSE) & ((match.template %in% rownames(data)) == FALSE)){
        stop('when finding a match for genes, match.template must be set to either a gene name (in the rownames of input data) or
             be the name of an annotation found in the annotations dataframe in the params list object, if neither is applicable,
             please use custom.template')
      }

      if (match.template %in% colnames(temp.annotations)) {
        if (sum(!is.na(temp.annotations)) != 0) {
          if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
            stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
          }
          temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]}

        template <- rep(0,nrow(temp.annotations))
        template[which(temp.annotations[,match.template] %in% annotation.level.set.high)] <- 1
      }
      if (match.template %in% rownames(data)) {
        template <- data[which(rownames(data) == match.template),]
      }

    }else{ template <- custom.template}

    p.vals <- apply(data,1,function(x)(cor.test(template,x,method=method, use=NA.handling)$p.value))
    corrs <- apply(data,1,function(x)(cor.test(template,x,method=method, use=NA.handling)$estimate))

    if(cut.by=="pvals"){corrs.past.cutoff <- names(p.vals)[which(p.vals < cutoff & corrs > 0)]}
    if(cut.by=="rvals"){corrs.past.cutoff <- names(corrs)[which(corrs > cutoff)]}
  }

  if (tolower(Find.Match.For) == "samples"){

    temp.annotations.genes <- params$annotations.genes



    if(is.na(custom.template)==TRUE){

      if (((match.template %in% colnames(temp.annotations)) == FALSE) & ((match.template %in% rownames(data)) == FALSE)){
        stop('when finding a match for samples, match.template must be set to either a sample name (in the colnames of input data) or
             be the name of an annotation found in the annotations.genes dataframe in the params list object, if neither is applicable,
             please use custom.template')
      }



      if (match.template %in% colnames(temp.annotations.genes)) {

        if (sum(!is.na(temp.annotations.genes)) != 0) {
          if (sum(rownames(data) %notin% rownames(temp.annotations.genes)) != 0 ) {
            stop('rownames of input data do not match rownames of annotations, cannot link annotations to data')
          }
          temp.annotations.genes <- temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),, drop = FALSE]}


        template <- rep(0,nrow(temp.annotations.genes))
        template[which(temp.annotations.genes[,match.template] %in% annotation.level.set.high)] <- 1
      }
      if (match.template %in% colnames(data)) {
        template <- data[,which(colnames(data) == match.template)]
      }

    }else{ template <- custom.template}

    p.vals <- apply(data,2,function(x)(cor.test(template,x,method=method, use=NA.handling)$p.value))
    corrs <- apply(data,2,function(x)(cor.test(template,x,method=method, use=NA.handling)$estimate))


    if(cut.by=="pvals"){corrs.past.cutoff <- names(p.vals)[which(p.vals < cutoff & corrs > 0)]}
    if(cut.by=="rvals"){corrs.past.cutoff <- names(corrs)[which(corrs > cutoff)]}
  }
  if(return.vals==FALSE){return(corrs.past.cutoff)}
  if(return.vals==TRUE){cbind(p.vals,corrs)}
}




#' @export
corrs2Gene <- function(  ##find correlations to specific gene, if no limits are provided, will show a histogram, with limits with create heatmap of genes passing correlation cutoff
  data, ##same data as before, should have genes in rows
  gene, ##gene name of gene you want correlations to
  limits =NULL, ## vector of two numbers, giving the lower and upper bounds
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  linkage = "complete",   ##linkage method, defulat complete linkage, can be changed
  NA.handling = "pairwise.complete.obs",   ##use for correlations, can be overwritten
  nbreaks=20,
  show.report=F,
  clust.rows = T, ##default for clustering rows, can be overwritten if error
  clust.cols = T, ##same as clust.cols but for cols
  row.groups = NA, ##number of groups to break the rows into based on dendrogram, can be overwritten
  col.groups = NA, ##same as col.groups but for cols, can be overwritten
  gaps.row = NULL, ##list of where to cut the rows if they're not clustered
  gaps.col = NULL, ##same as gaps.row but for columns
  gap.width=1,
  order.by.gene = NULL, ##gene to order by
  order.by.sample = NULL, ##sample to order by
  cell.width = NA,
  cell.height = NA,
  fontsize.row = 10,
  fontsize.col = 10,
  show.rownames=T,
  show.colnames=F,
  treeheight.row=20,
  treeheight.col=20,
  hide.plot=FALSE,
  na.fix=FALSE,
  na.offset =2,
  show.legend=TRUE,
  show.annotations=TRUE,
  is.raw.Ct=FALSE, ##if true, will reverse color scale to show yellow as high expressing
  drop.annot.levels = TRUE
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  main <- paste("Correlations_to_",gene)
  gen.cor<-cor(x=(data[which(rownames(data)==gene),]),y=t(data), use=NA.handling, method=method)

  if(is.null(limits)==TRUE){
    hist(gen.cor,breaks=nbreaks, main=main)
  }

  if(is.null(limits)==FALSE){

    below <- limits[1]
    above <- limits[2]

    subset<-data[which(gen.cor>above| gen.cor< below),];


    myHeatmap(subset, method = method, linkage = linkage, NA.handling = NA.handling, clust.rows = clust.rows, clust.cols = clust.cols, row.groups = row.groups,
              gaps.row = gaps.row, gaps.col = gaps.col, gap.width= gap.width, order.by.gene = order.by.gene, order.by.sample = order.by.sample,
              cell.width = cell.width, cell.height = cell.height, fontsize.row = fontsize.row, fontsize.col = fontsize.col, show.rownames = show.rownames, show.colnames = show.colnames,
              treeheight.row = treeheight.row, treeheight.col = treeheight.col, hide.plot = hide.plot, na.fix = na.fix,na.offset = na.offset, show.legend = show.legend,
              show.annotations = show.annotations, is.raw.Ct = is.raw.Ct, drop.annot.levels = drop.annot.levels, main = paste(main, "Below:", below,"Above:", above))

    if(show.report==T){ report <- gen.cor[,which(gen.cor>above| gen.cor< below)]; return(report)}

  }

}





#' @export
correlateGenes <- function(  ##broad gene correlations
  data, ##data matrix with genes in the rows
  limits =NULL, ## vector of two numbers, giving the lower and upper bounds
  nbreaks=20,
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  NA.handling = "pairwise.complete.obs"   ##use for correlations, can be overwritten
){
  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  main=paste("Gene Correlations")
  cor.dat <- cor(t(data),use=NA.handling,method=method)
  cor.dat[!upper.tri(cor.dat)] <- NA
  if(is.null(limits)==TRUE){
    hist(cor.dat,breaks=nbreaks, main=main)
    #return(cor.dat)
  }

  if(is.null(limits)==FALSE){

    below <- limits[1]
    above <- limits[2]


    picked.cors <- which((cor.dat > above & cor.dat < .99)| cor.dat < below ,arr.ind = T)

    cor.genes <- as.matrix(apply(picked.cors,2,function(x)(colnames(t(data))[x])))
    genes.and.cors <- as.data.frame(cbind(cor.genes, Correlation=cor.dat[picked.cors]))
    colnames(genes.and.cors) <- c("Gene1","Gene2","Correlation")
    return(genes.and.cors)
  }
}



#' @export
reportGenes <- function(  ##returns report summary of range of expression
  data, ##dataset, genes should be in rows
  list, ##list of genes to summarize
  exact = T,
  ranges="fixed", ##default, can be changed, other option is mad
  fixed.range =2, #only if ranges=fixed, can be changed to adjust range
  weight = 1.25 #only used if ranges=mad, weight to be added to mad for ranges
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  ##subset for list
  if (exact == TRUE) {
    list <- list
    if (sum(list %in% rownames(data)) == 0 ) {stop('exact matches for list not found in rownames data')}
  }else{
    list <- rownames(data)[grep(paste(list, collapse = "|"),rownames(data))]
    if (sum(list %in% rownames(data)) == 0 ) {stop('inexact matches for list not found in rownames data')}
  }

  report <- matrix(nrow=length(list),ncol=6)

  for (i in 1:length(list)) {

    gene <- data[which(rownames(data) %in% list[i]),]
    percentNA <- sum(is.na(gene))/length(gene)
    percentdetected <- sum(!is.na(gene))/length(gene)
    gene <- gene[which(!is.na(gene))]
    median <- median(gene, na.rm=T)

    if (ranges=="fixed") {
      upper <- median - fixed.range
      lower <- median + fixed.range
    }

    if (ranges == "mad") {
      upper <- median - (mad(gene,na.rm=T)*weight)
      lower <- median + (mad(gene,na.rm=T)*weight)
    }

    high <- round(sum(gene < upper)/length(gene),2)
    middle <- round(sum(gene > upper & gene < lower)/length(gene),2)
    low <- round(sum(gene > lower)/length(gene),2)

    report[i,1] <- list[i]
    report[i,2] <- percentdetected
    report[i,3] <- percentNA
    report[i,4] <- high
    report[i,5] <- middle
    report[i,6] <- low
  }

  colnames(report) <- c("Gene","Percent_Samples_Detected","Percent_Samples_NA","Percent_High_Expressing",
                        "Percent_Mid-Range_Expression","Percent_Low_Expressing")

  return(as.data.frame(report))
}



#' @import ggplot2
#' @importFrom cluster silhouette
#' @importFrom gtools permute
#' @export
find.silhouette <- function(
  data,
  to.sil = "samples", # can change to genes
  to.view="rand.to.clust", ##can also be "all.clusts" or "rand.all.clusts"
  ngroups, #for rand to clust
  maxgroups=12,
  max.iter=10,
  method = "pearson",
  NA.handling = "pairwise.complete.obs",
  linkage = "complete",
  main="Average Silhouette Width",
  axis.label="Silhouette Width", #x axis for rand.to clust, y label for the others
  main.label.size= 30,
  axis.label.size=20,
  legend.position = "bottom"
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  ###clustering
  if (method %in% c("spearman","pearson", "kendall")) {
    clust.genes<-(as.dist(1-cor(t(data),method=method,use=NA.handling)));
    clust.samps<-(as.dist(1-cor(data,method=method,use=NA.handling)))
  }


  if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
    clust.genes <- dist(data,method=method)
    clust.samps <- dist(t(data), method=method)
  }


  ##cluster samples and cut
  if (to.sil=="samples") {clust.mat <- clust.samps}
  if (to.sil=="genes") {clust.mat <- clust.genes}

  clusts <- hclust(clust.mat, method = linkage)

  if (to.view =="rand.to.clust") {   ##look at specific number of clusters to randomized

    groupings <- cutree(clusts, k=ngroups)

    sil.coef <- summary(cluster::silhouette(groupings, clust.mat))$avg.width

    rands <-  replicate(max.iter,summary(cluster::silhouette(gtools::permute(groupings), clust.mat))$avg.width, simplify = "vector")

    data.to.plot <- data.frame(vals=c(rands, sil.coef))

    p <- ggplot(data.to.plot, aes(x=vals)) + geom_density(size=1) + geom_point(x=sil.coef, y=density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x-sil.coef))], size=3) + geom_text(x=sil.coef, y=(density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x-sil.coef))])+5, label=round(sil.coef,3)) +
      geom_segment(x=sil.coef, xend=sil.coef, y=0, yend=density(data.to.plot$vals)$y[which.min(abs(density(data.to.plot$vals)$x-sil.coef))], size=1) + theme_bw() +  theme(panel.grid=element_blank() ,
                                                                                                                                                                           axis.title = element_text(size=axis.label.size), plot.title = element_text(size=main.label.size, hjust = 0.5))+
      xlab(paste(axis.label)) + ylab("Density") + ggtitle(paste(main))


  }

  ###if checking for how many clusters

  if (to.view != "rand.to.clust") {


    sil.coefs <- c()


    for (nclust in 2:maxgroups) {
      groupings <- cutree(clusts, k=nclust)
      sil.coefs <- c(sil.coefs, summary(cluster::silhouette(groupings, clust.mat))$avg.width)

    }

    if (to.view == "all.clusts") {
      to.plot <- data.frame(clusts=2:maxgroups, coefs=sil.coefs, type="Original Data")
      p <- ggplot(to.plot,aes(x=clusts, y=sil.coefs)) + geom_line(size=1) + geom_point(size=3) + theme_bw()+ theme(panel.grid = element_blank()) + xlab("Clusters") +
        ylab(paste(axis.label)) + ggtitle(paste(main)) + theme(axis.title = element_text(size=axis.label.size), plot.title = element_text(size=main.label.size, hjust=0.5))+ scale_x_continuous(breaks=scales::pretty_breaks())
    }

    if (to.view == "rand.all.clusts") {

      rands <- NULL

      for (nclust in 2:maxgroups) {
        groupings <- cutree(clusts, k=nclust)
        rands <- cbind(rands,replicate(max.iter,summary(cluster::silhouette(gtools::permute(groupings), clust.mat))$avg.width, simplify = "vector"))

      }


      to.plot <- data.frame(clusts=2:maxgroups, coefs=sil.coefs, type="Original Data")
      colnames(rands) <- 2:maxgroups
      rands.melt <- reshape2::melt(rands); rands.melt$type <- "Randomized Data"



      p <- ggplot() + geom_line(data=to.plot,aes(x=clusts, y=sil.coefs),size=1) + geom_point(data=to.plot,aes(x=clusts, y=sil.coefs, shape=type),size=3) +
        geom_boxplot(data=rands.melt, aes(x=Var2, y=value, group=Var2)) + geom_point(data=rands.melt, aes(x=Var2,y=value, group=Var2, shape=type), size=1) +
        stat_summary(data=rands.melt,aes(x=Var2, y=value), fun = median, geom = 'line') +
        stat_summary(data=rands.melt,aes(x=Var2, y=value), fun = median, geom = 'point', size=4,shape=17) +
        theme_bw()+ theme(panel.grid = element_blank()) + xlab("Clusters") + scale_shape_manual(values=c(19,17),labels=c("Original Data", "Randomized Data")) +
        ylab(paste(axis.label)) + ggtitle(paste(main)) + theme(axis.title = element_text(size=axis.label.size), plot.title = element_text(size=main.label.size, hjust=0.5),
                                                               legend.position = legend.position, legend.direction="horizontal",legend.title = element_blank(), legend.background = element_rect(color="black")) + scale_x_continuous(breaks=scales::pretty_breaks()) +
        ylim(c(min(rands), max(sil.coefs)))

    }
  }
  return(p)
}

