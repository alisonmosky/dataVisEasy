
###Initialize and Update Parameters####
parameters <- list(scale.range = c(-1,1),
               scale.colors = c("blue","black","yellow"),
               n.colors.range = 100,
               annotations = NA,
               annot_samps = NA,
               annotations.genes = NA,
               annot_genes = NA,
               annot_cols = NA,
               expression_gradient.colors = c("blue","lightblue","gray","indianred","firebrick"))


initiate_params <- function(params){
  assign("params",params,envir =  .GlobalEnv)
}


set_scale.range <- function(range){
  if(length(range) == 2 & class(range) == "numeric" & range[2] > range[1]){
    params$scale.range <<- range}
  else{stop('range must be a numeric of length 2 and the second element must be greater than the first')}}

set_scale.colors <- function(colors){#unlockBinding("params", env = as.environment("package:dataVisEasy"));
  params$scale.colors <<- colors}

set_n.colors.range <- function(n){#unlockBinding("params", env = as.environment("package:dataVisEasy"));
  params$n.colors.range <<- n}


set_annotations <- function(annotations){
  params$annotations <<- annotations
  if (sum(is.na(params$annotations) != 0)) {
    if (length(annotations) ==1 ) {
      params$annotations <<- NA
    }else{
    params$annotations[is.na(params$annotations)] <<- "No_Annot"
    warning('NAs present in annotations, converted to "No_Annot" for functionality purposes')
    }
  }
}

update_annotations <- function(   ###should make an option to full join if stuff thats being added doesn't cover everything....
  annotation,
  values
){
  if (annotation %in% colnames(annotations)) {
    params$annotations[,annotation] <<- values
  }else{
    params$annotations$V1 <<- values
    colnames(params$annotations)[colnames(params$annotations)=="V1"] <<- annotation
  }
}

set_annot_samps <- function(annotations= NULL){
  if (is.null(annotations) == TRUE) { params$annot_samps <<- params$annotations
  }else{
    suppressWarnings(if (is.na(annotations)) { params$annot_samps <<- NA
    }else{
        params$annot_samps <<- params$annotations[,which(colnames(params$annotations) %in% annotations), drop = FALSE]
        params$annot_samps <<- params$annot_samps[,match(annotations, colnames(params$annot_samps)), drop = FALSE]
    })
  }
}


set_annotations.genes <- function(annotations.genes){
  params$annotations.genes <<- annotations.genes
  if (sum(is.na(params$annotations.genes) != 0)) {
    if (length(annotations.genes) == 1) {
      params$annotations.genes <<- NA
    }else {
    params$annotations.genes[is.na(params$annotations.genes)] <<- "No_Annot"

    warning('NAs present in annotations.genes, converted to "No_Annot" for functionality purposes')
    }
  }
}

update_annotations.genes <- function(  ##right now can only do one at a time
  annotation,
  values
){
  if (annotation %in% colnames(annotations.genes)) {
    params$annotations.genes[,annotation] <<- values
  }else{
    params$annotations.genes$V1 <<- values
    colnames(params$annotations.genes)[colnames(params$annotations.genes)=="V1"] <<- annotation
  }
}

set_annot_genes <- function(annotations = NULL){
  if (is.null(annotations) == TRUE) { params$annot_genes <<- params$annotations.genes
  }else{
    suppressWarnings( if (is.na(annotations)) { (params$annot_genes <<- NA)
    }else{
        params$annot_genes <<- params$annotations.genes[,which(colnames(params$annotations.genes) %in% annotations),drop=FALSE]
        params$annot_genes <<- params$annot_genes[,match(annotations, colnames(params$annot_genes)), drop = FALSE]
    })

  }

}



set_annot_cols <- function(annot_cols){
  params$annot_cols <<- annot_cols
  suppressWarnings( if (is.na(annot_cols)) { params$annot_cols <<- NA})
}

update_annot_cols <- function(  ##right now must be one at a time
  annotation,
  values.list
){
  if (annotation %in% names(params$annot_cols)) {
    params$annot_cols[annotation] <<- values.list
  }else{
    params$annot_cols$V1 <<- values.list
    names(params$annot_cols)[names(params$annot_cols)=="V1"] <<- annotation
  }
}


set_expression_gradient.colors <- function(colors){
  if (length(colors) == 5) {
    is.color = sapply(colors, function(x) {
      tryCatch( is.matrix( col2rgb( x ) ), error = function( e ) FALSE )
    })
    if(any(is.na(names(is.color)))) is.color[is.na(names(is.color))] = FALSE

    if (sum(is.color) == 5) {
      params$expression_gradient.colors <<- colors
    }else{stop('must contain valid colors')}
  }else{stop('must be a vector of colors of length 5')}}

#########
"%notin%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)

assessScale <- function(     ###Returns percent of data above, below, and within range of scale for heatmap
  data ##data to be assessed,
  ##scale.range, range of heatmap
){
  above <- sum(data > params$scale.range[2], na.rm=T)/sum(!is.na(data))
  within <- sum(data <= params$scale.range[2] & data >= params$scale.range[1], na.rm=T)/sum(!is.na(data))
  below <- sum(data < params$scale.range[1], na.rm=T)/sum(!is.na(data))
  report <- c(below, within, above); names(report) <- c("Percent Below Range", "Percent Within Range","Percent Above Range")
  return(report)
}

subsetGenes <- function(   ##simply pulls list of genes out of dataset into a new matrix
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list,   ##list of gene or genes that will be extracted from the dataset, exact matches
  order.by = NULL, ##gene to order by
  exact = TRUE ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
){

  if(exact==TRUE){
    subset <- data[which(rownames(data) %in% list),];
    if (length(subset) == 0 ) {stop('exact matches for list not found in rownames data')}
  }

  if (exact==FALSE){
    subset <- data[grep(paste(list, collapse = "|"),rownames(data)),]
    if (length(subset) == 0 ) {stop('inexact matches for list not found in rownames data')}
  }

  if (is.null(order.by)==FALSE){
    subset <- subset[,order(data[which(rownames(data) %in% order.by),],na.last = F)]
  }
  return(subset)
}


##need to fix this to go with annotations
subsetSamples <- function(  ##subset samples out of matrix based on metadata group
  data, ##matrix of values, with genes in the rows and samples in the columns
  group, ##group from which to subset, vector that is same length as columns
  take.out ##which part of the group to extract, one or more of the factors in the "group"
){

  temp.annotations <- params$annotations

  if (group %in% colnames(temp.annotations)) {

    if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
      stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
    }
    temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]

    groupings <- as.factor(temp.annotations[,group] )

    if (sum(take.out %in% groupings) != length(take.out)) {stop('provided arguments to take.out not found in indicated group')}
    set <- data[,which(groupings %in% take.out)]
  }else{
    if (sum(take.out %in% group) != length(take.out)) {stop('provided arguments to take.out not found in indicated group')}
    set <- data[,which(group %in% take.out)]
  }
  return(set)
}

myColorRamp5 <- function(colors, values, percent.mad = 0.5) {  ###color data over a range, assumes 5 colors, sets in quadrants according to median +/- mad
  out <- rep(rgb(0,0,0),length(values))
  for(i in 1:length(values)){
    if(is.na(values[i])){
    } else{
      if (values[i] <= (median(values,na.rm=T)-percent.mad*mad(values,na.rm=T))){
        v <- (values[i] - min(values,na.rm=T))/( (median(values,na.rm=T)-percent.mad*mad(values,na.rm=T)) - min(values, na.rm=T) )
        x <- colorRamp(colors[1:2])(v)
        out[i] <- rgb(x[,1], x[,2], x[,3], maxColorValue = 255)}

      if (values[i]>median(values,na.rm=T)-percent.mad*mad(values,na.rm=T) & values[i]<=median(values, na.rm=T)){
        v <- (values[i] - (median(values,na.rm=T)-percent.mad*mad(values,na.rm=T)) )/ (median(values, na.rm=T) -(median(values,na.rm=T)-percent.mad*mad(values,na.rm=T)))
        x <- colorRamp(colors[2:3])(v)
        out[i] <- rgb(x[,1], x[,2], x[,3], maxColorValue = 255)}

      if (values[i]<=median(values,na.rm=T)+percent.mad*mad(values,na.rm=T) & values[i]>median(values, na.rm=T)){
        v <- (values[i] - median(values, na.rm=T))/ ( (median(values,na.rm=T)+percent.mad*mad(values,na.rm=T))- median(values, na.rm=T))
        x <- colorRamp(colors[3:4])(v)
        out[i] <- rgb(x[,1], x[,2], x[,3], maxColorValue = 255)}

      if (values[i]>median(values,na.rm=T)+percent.mad*mad(values,na.rm=T)){
        v <- (values[i] - (median(values,na.rm=T)+percent.mad*mad(values,na.rm=T)))/(max(values, na.rm=T)  - (median(values,na.rm=T)+percent.mad*mad(values,na.rm=T)))
        x <- colorRamp(colors[4:5])(v)
        out[i] <- rgb(x[,1], x[,2], x[,3], maxColorValue = 255)}
    }
  }
  return(out)
}


myHeatmap <- function(  ##basic heatmap, can subset for gene list
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list = NULL,   ##list of gene or genes that will be extracted from the dataset, exact matches
  exact = TRUE, ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  linkage = "complete",   ##linkage method, defulat complete linkage, can be changed
  NA.handling = "pairwise.complete.obs",   ##use for correlations, can be overwritten
  clust.rows = T, ##default for clustering rows, can be overwritten if error
  clust.cols = T, ##same as clust.cols but for cols
  row.groups = NA, ##number of groups to break the rows into based on dendrogram, can be overwritten
  col.groups = NA, ##same as col.groups but for cols, can be overwritten
  gaps.row = NULL, ##list of where to cut the rows if they're not clustered
  gaps.col = NULL, ##same as gaps.row but for columns
  gap.width=1,
  main = NULL,  ##for title of chart, must be in quotes, will default to list of the genes
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
  na.offset = 2,
  show.legend=TRUE,
  show.annotations=TRUE,
  is.raw.Ct=FALSE, ##if true, will reverse color scale to show yellow as high expressing
  drop.annot.levels=TRUE
){

  if (is.null(main)==TRUE){
    main <- paste("Genes of Interest:",paste(list, collapse = ","))}

  if (is.null(list) == TRUE) {list <- rownames(data)}


  ##subset for list
  if (exact == TRUE) {
    subset <- data[which(rownames(data) %in% list),]
    if (length(subset) == 0 ) {stop('exact matches for list not found in rownames data')}
  }else{
    subset <- data[grep(paste(list, collapse = "|"),rownames(data)),]
    if (length(subset) == 0 ) {stop('inexact matches for list not found in rownames data')}
  }


  if (na.fix==TRUE) {
    if (is.raw.Ct==TRUE) {subset[which(is.na(subset))] <- max(subset,na.rm=T)+na.offset}
    if (is.raw.Ct==FALSE) {subset[which(is.na(subset))] <- min(subset,na.rm=T)-na.offset}
  }


  if (method %in% c("spearman","pearson", "kendall")) {
    clust.genes<-(as.dist(1-cor(t(subset),method=method,use=NA.handling)));
    clust.samps<-(as.dist(1-cor(subset,method=method,use=NA.handling)))
  }

  if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
    clust.genes <- dist(subset,method=method)
    clust.samps <- dist(t(subset), method=method)
  }


  if (is.null(order.by.gene)==FALSE){
    if(is.raw.Ct==FALSE){subset <- subset[,order(data[which(rownames(data) %in% order.by.gene),],na.last = F)]}
    if(is.raw.Ct==TRUE){subset <- subset[,order(data[which(rownames(data) %in% order.by.gene),],na.last = T)]}
    clust.cols <- F
  }

  if (is.null(order.by.sample)==FALSE){
    if(is.raw.Ct==FALSE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = F),]}
    if(is.raw.Ct==TRUE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = T),]}
    clust.rows <- F
  }

  if(clust.rows==T){heightrow <- treeheight.row}
  if(clust.cols==T){heightcol <- treeheight.col}

  subset1 <- subset
  subset <- scales::squish(subset,params$scale.range)
  breaks <- seq(params$scale.range[1], params$scale.range[2],length.out=params$n.colors.range)
  my_cols=colorRampPalette(params$scale.colors)(n=params$n.colors.range-1)
  if(is.raw.Ct==TRUE){my_cols <- rev(my_cols)}

  if(na.fix==TRUE){
    if(is.raw.Ct==TRUE){
      subset[which(subset1==max(subset1))] <- params$scale.range[2]+0.04
      breaks <- c(breaks, params$scale.range[2]+0.01, params$scale.range[2]+0.05)
      my_cols <- c(my_cols,params$scale.colors[1],"grey90")      ##may need an option to set na_col
    }
    if(is.raw.Ct==FALSE){
      subset[which(subset1==min(subset1))] <- params$scale.range[1]-0.04
      breaks <- c(params$scale.range[1]-0.05,params$scale.range[1]-0.01,breaks)
      my_cols <- c("grey90",params$scale.colors[1],my_cols)
    }
  }


  temp.annot_samps <- params$annot_samps
  temp.annot_genes <- params$annot_genes
  temp.annot_cols <- params$annot_cols

  if (drop.annot.levels == TRUE) {
    suppressWarnings( if (is.na(temp.annot_samps) == F) {
      temp.annot_samps[] <- lapply(temp.annot_samps, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_samps <- temp.annot_samps %>% tibble::rownames_to_column("Sample")
      temp.annot_samps <- droplevels(temp.annot_samps[which(temp.annot_samps$Sample %in% colnames(subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="Sample")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_samps) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0 ) {
        for (annot.i in 1:length(spec.cols)) {
          annot <- colnames(temp.annot_samps)[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_samps[,which(colnames(temp.annot_samps)==annot)])  )]
          if ( sum( levels(temp.annot_samps[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_samps[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_samps[,annot])[levels(temp.annot_samps[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
            }
          }
      }
    })


    suppressWarnings( if (is.na(params$annot_genes) == F) {
      temp.annot_genes[] <- lapply(temp.annot_genes, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_genes <- temp.annot_genes %>% tibble::rownames_to_column("Gene")
      temp.annot_genes <- droplevels(temp.annot_genes[which(temp.annot_genes$Gene %in% rownames(subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="Gene")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_genes) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0) {
        for (annot.i in 1:length(colnames(temp.annot_genes))) {
          annot <- colnames(temp.annot_genes)[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_genes[,which(colnames(temp.annot_genes)==annot)])  )]
          if ( sum( levels(temp.annot_genes[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_genes[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_genes[,annot])[levels(temp.annot_genes[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    })

  }




  if (clust.cols == T) {
    tryclustcols <- try(hclust(clust.samps, linkage), silent = T)
    if (class(tryclustcols) == "try-error") {stop('cannot cluster columns, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.cols = F or specify order.by.gene')}
    }

  if (clust.rows == T) {
    tryclustrows <- try(hclust(clust.genes, linkage), silent = T)
    if (class(tryclustrows) == "try-error") {stop('cannot cluster rows, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.rows = F or specify order.by.sample')}
    }



  pheatmap(subset,col=my_cols, breaks=breaks, border_color = NA, clustering_method=linkage,annotation_col=temp.annot_samps, annotation_colors = temp.annot_cols,
           clustering_distance_rows = clust.genes, clustering_distance_cols = clust.samps, main=paste(main,"\n Method_",method,"_Linkage_",linkage),
           cluster_rows = clust.rows, cluster_cols = clust.cols, cutree_rows = row.groups, cutree_cols = col.groups, gaps_row = gaps.row, gaps_col = gaps.col,
           cellwidth = cell.width, cellheight = cell.height, fontsize_row = fontsize.row, fontsize_col = fontsize.col, show_rownames = show.rownames,show_colnames = show.colnames,
           treeheight_row = heightrow ,treeheight_col = heightcol, silent = hide.plot, legend=show.legend, annotation_legend = show.annotations,
           annotation_row=temp.annot_genes, drop_levels = drop.annot.levels)


    }



myHeatmapByAnnotation <- function(
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list = NULL,   ##list of gene or genes that will be extracted from the dataset, exact matches
  exact = TRUE, ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
  groupings, ##either character vector pointing to annotations, dataframe where the first row will be taken, or factor, if unnamed factor, wont properly order and will assume in same order as data
  groupings.gaps = NULL,
  groupings.genes = FALSE,
  groupings.genes.gaps = NULL,
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  linkage = "complete",   ##linkage method, defulat complete linkage, can be changed
  NA.handling = "pairwise.complete.obs",   ##use for correlations, can be overwritten
  clust.rows = TRUE, ##default for clustering rows, can be overwritten if error
  clust.cols = TRUE, ##same as clust.cols but for cols
  row.groups = NA, ##number of groups to break the rows into based on dendrogram, can be overwritten
  col.groups = NA, ##same as col.groups but for cols, can be overwritten
  gaps.row = TRUE, ##list of where to cut the rows if they're not clustered
  gaps.row.spec = NULL,
  gaps.col = TRUE, ##same as gaps.row but for columns
  gaps.col.spec = NULL, ##Null if want to separate automatically by group, can override with vector of indices to split
  gap.width=1, ##width of gaps in between
  main = NULL,  ##for title of chart, must be in quotes, will default to list of the genes
  order.by.gene = NULL, ##gene to order by
  order.by.sample=NULL,
  fontsize.row = 10,
  fontsize.col = 10,
  show.rownames=T,
  show.colnames=F,
  treeheight.row=20,
  treeheight.col=20,
  cell.width=NA, ##can change cell width
  cell.height=NA, ##can change cell height
  hide.plot=FALSE,
  na.fix=FALSE,
  na.offset = 2,
  is.raw.Ct=FALSE,
  show.legend=TRUE,
  show.annotations=TRUE,
  drop.annot.levels = TRUE
){


  if (is.null(main)==TRUE){
    main <- paste("Genes of Interest:",paste(list, collapse = ","))}

  if (is.null(list) == TRUE) {list <- rownames(data)}

  ##subset for list
  if (exact == TRUE) {
    data.subset <- as.matrix(data[which(rownames(data) %in% list),]);colnames(data.subset) <- colnames(data)
    if (length(data.subset) == 0 ) {stop('exact matches for list not found in rownames data')}
    #if (groupings.genes[which(!is.na(groupings.genes))[1]] != FALSE) {groupings.genes <- droplevels(groupings.genes[which(rownames(data) %in% list)])}

  }else{
    data.subset <- as.matrix(data[grep(paste(list, collapse = "|"),rownames(data)),]); colnames(data.subset) <- colnames(data)
    if (length(data.subset) == 0 ) {stop('exact matches for list not found in rownames data')}
    #if (groupings.genes[which(!is.na(groupings.genes))[1]] != FALSE ) {groupings.genes <- droplevels(groupings.genes[grep(paste(list, collapse = "|"), rownames(data))]) }
  }


  ####new code to order annotations if not in order or if extra etc####
  suppressWarnings( if (is.na(params$annotations) == FALSE) {
    if (sum(colnames(data.subset) %notin% rownames(params$annotations)) != 0 ) {
      stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
    }
    temp.annotations <- params$annotations[match(colnames(data.subset), rownames(params$annotations)),, drop = FALSE]}  )
  suppressWarnings( if (is.na(params$annotations.genes) == FALSE) {
    if (sum(rownames(data.subset) %notin% rownames(params$annotations.genes)) != 0 ) {
      stop('rownames of input data do not match rownames of annotations.genes, cannot link annotations to data')
    }
    temp.annotations.genes <- params$annotations.genes[match(rownames(data.subset), rownames(params$annotations.genes)),, drop = FALSE]}  )




  if (class(groupings) =="character"){
    if ( sum(groupings %notin% colnames(temp.annotations)) != 0) {
      stop('supplied character vector for groupings not found in sample annotations')}
    factorgroupings <- makefactorgroup(temp.annotations, groupings, specify.gaps = groupings.gaps, return.gaps = TRUE)
    groupings <- factorgroupings$factor.group
    gaps.groupings <- c(factorgroupings$gaps)
  }else{
    if (class(groupings) == "data.frame") {
      ##order groupings by order of subset
      groupings <- droplevels(groupings[match(colnames(data.subset), rownames(groupings)),1]) %>% as.factor()
      if (sum(is.na(groupings) != 0 )) {
        groupings <- as.character(groupings); groupings[is.na(groupings)] <- "No_Annot"
        groupings <- as.factor(groupings)}
    }
    if (class(groupings) == "factor"){
      if (is.null(names(groupings))==FALSE) {
        groupings <- droplevels(groupings[match(colnames(data.subset), names(groupings))]) %>% as.factor()
        if (sum(is.na(groupings) != 0 )) {
          groupings <- as.character(groupings); groupings[is.na(groupings)] <- "No_Annot"
          groupings <- as.factor(groupings)}
      }else{
        if (length(groupings) != ncol(data.subset)) {
          stop('unnamed factor supplied to groupings not the same length as number of columns in data')
        }
      }  ##do nothing and leave as is
    }
  }



  if (class(groupings.genes) =="character") {
    if (sum(groupings.genes %notin% colnames(temp.annotations.genes)) != 0) {
      stop('supplied character vector for groupings.genes not found in gene annotations')}
    factorgroupings.genes <- makefactorgroup(temp.annotations.genes, groupings.genes, specify.gaps = groupings.genes.gaps, return.gaps = TRUE)
    groupings.genes <- factorgroupings.genes$factor.group
    gaps.groupings.genes <- c(factorgroupings.genes$gaps)
  }else{
    if (class(groupings.genes) == "data.frame") {
      ##order groupings by order of subset
      groupings.genes <- droplevels(groupings.genes[match(rownames(data.subset), rownames(groupings.genes)),1]) %>% as.factor()
      if (sum(is.na(groupings.genes) != 0 )) {
        groupings.genes <- as.character(groupings.genes); groupings.genes[is.na(groupings.genes)] <- "No_Annot"
        groupings.genes <- as.factor(groupings.genes)}
    }
    if (class(groupings.genes) == "factor"){
      if (is.null(names(groupings.genes))==FALSE) {
        groupings.genes <- droplevels(groupings.genes[match(rownames(data.subset), names(groupings.genes))]) %>% as.factor()
        if (sum(is.na(groupings.genes) != 0 )) {
          groupings.genes <- as.character(groupings.genes); groupings.genes[is.na(groupings.genes)] <- "No_Annot"
          groupings.genes <- as.factor(groupings.genes)}
      }else{
        if (length(groupings.genes) != nrow(data.subset)) {
          stop('unnamed factor supplied to groupings.genes not the same length as number of rows from data or number of rows after subsetting for supplied list')
        }
      }  ##do nothing and leave as is
    }
  }





  ####end of new code to order groupings




  if (na.fix==TRUE) {
    if (is.raw.Ct==TRUE) {data.subset[which(is.na(data.subset))] <- max(data,na.rm=T)+na.offset}
    if (is.raw.Ct==FALSE) {data.subset[which(is.na(data.subset))] <- min(data,na.rm=T)-na.offset}
  }

  if (is.null(order.by.gene) == FALSE) {if ((order.by.gene %in% rownames(data.subset)) == FALSE) { order.by.gene <- NULL}}

  if (is.null(order.by.sample) == FALSE) {if ((order.by.sample %in% colnames(data.subset)) == FALSE) { order.by.sample <- NULL}}



  ###groupings


  if (sum(groupings != FALSE, na.rm=T) != 0) {

    ind.col=0     ##needed if not connected to annotations
    gaps.col.n = NULL

    samp.order <- NULL

    for (i in 1:length(unique(groupings))) {
      subset <- as.matrix(data.subset[,which(groupings==levels(groupings)[i])]); colnames(subset) <- colnames(data.subset)[which(groupings==levels(groupings)[i])] #as.matrix is necessary if there is a group of one


      ##cluster samples of each group
      if (ncol(subset) > 1) {
        if(clust.cols==TRUE){

          if (method %in% c("spearman","pearson", "kendall")) {
            clust.samps<-(as.dist(1-cor(subset,method=method,use=NA.handling)))
          }

          if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
            clust.samps <- dist(t(subset), method=method)
          }


          tryclustcols <- try(hclust(clust.samps, linkage), silent = T)
          if (class(tryclustcols) == "try-error") {stop('cannot cluster columns, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.cols = F or specify order.by.gene')}

          a <- hclust(clust.samps, method = linkage)
          samp.order <- c(samp.order,a$labels[a$order])}
        if(clust.cols==FALSE){
          samp.order <- c(samp.order,colnames(subset))
        }
        }
      if(ncol(subset)==1){samp.order <- c(samp.order,colnames(data.subset)[which(groupings==levels(groupings)[i])])}



      ##order genes for each sample group if necessary
      if (is.null(order.by.gene)==FALSE){
        if(is.raw.Ct==FALSE){subset <- subset[,order(subset[which(rownames(subset) %in% order.by.gene),],na.last = F)]}
        if(is.raw.Ct==TRUE){subset <- subset[,order(subset[which(rownames(subset) %in% order.by.gene),],na.last = T)]}
        clust.cols <- F

      }


      if (i==1) {combined <- subset}
      if (i!=1) {combined <- cbind(combined,subset)}


      if (exists("factorgroupings") == FALSE) {
        if(gaps.col==TRUE){
          ind.col <- ind.col + ncol(subset)
          gaps.col.n <- c(gaps.col.n,ind.col)
        }

      }

      }

    if (is.null(order.by.gene)==TRUE){combined <- combined[,samp.order]} #this is setting it, the order is determined in the loop, kind of redundant but dont feel like changing now

    if (sum(groupings.genes == FALSE, na.rm = T) != 0) {  ##if not going on to group genes, see if should be ordered by sample, set clustering of genes
      if (is.null(order.by.sample)==FALSE){
        if(is.raw.Ct==FALSE){combined <- combined[order(combined[,which(colnames(combined) %in% order.by.sample)],na.last = F),]}
        if(is.raw.Ct==TRUE){combined <- combined[order(combined[,which(colnames(combined) %in% order.by.sample)],na.last = T),]}
        clust.rows <- F
      }else {
        #clustering of genes
        if (method %in% c("spearman","pearson", "kendall")) {
          clust.genes<-(as.dist(1-cor(t(combined),method=method,use=NA.handling)))
        }
        if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
          clust.genes <- dist((combined), method=method)
        }
      }

      gaps.row <- NULL
    }




    if(gaps.col==TRUE){
      if (exists("factorgroupings") == FALSE) {
        gaps.col <- gaps.col.n[-length(gaps.col.n)]}else{gaps.col <- gaps.groupings}
      gaps.col <- sort(rep(gaps.col,gap.width))}

    clust.cols <- F

    data.subset <- combined ## if going on to gene groupings will start in the same format
  }


  ###replicate above and switch for rows




  if (sum(groupings.genes != FALSE, na.rm=T) != 0) {

    ind.row=0      ###needed if not pointing to annotations
    gaps.row.n = NULL

    gene.order <- NULL

    for (i in 1:length(unique(groupings.genes))) {
      subset <- as.matrix(data.subset[which(groupings.genes==levels(groupings.genes)[i]),]); #as.matrix is necessary if there is a group of one


      ##cluster samples of each group
      if (ncol(subset) > 1) {
        rownames(subset) <- rownames(data.subset)[which(groupings.genes==levels(groupings.genes)[i])]
        if(clust.rows==TRUE){

          #clustering of genes
          if (method %in% c("spearman","pearson", "kendall")) {
            clust.genes<-(as.dist(1-cor(t(subset),method=method,use=NA.handling)))
          }

          if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
            clust.genes <- dist((subset), method=method)
          }

          tryclustrows <- try(hclust(clust.genes, linkage), silent = T)
          if (class(tryclustrows) == "try-error") {stop('cannot cluster rows, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.rows = F or specify order.by.sample')}

          a <- hclust(clust.genes, linkage)
          gene.order <- c(gene.order,a$labels[a$order])}
        if(clust.rows==FALSE){
          gene.order <- c(gene.order,rownames(subset))
        }
        }
      if(ncol(subset)==1){
        colnames(subset) <- rownames(data.subset)[which(groupings.genes==levels(groupings.genes)[i])]
        gene.order <- c(gene.order,colnames(subset))
        subset <- t(subset)}



      #order samples for each gene group if necessary
      if (is.null(order.by.sample)==FALSE){
        if(nrow(subset) > 1) {
          if(is.raw.Ct==FALSE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = F),]}
          if(is.raw.Ct==TRUE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = T),]}
          clust.rows <- F
        }

      }


      if (i==1) {combined <- subset}
      if (i!=1) {combined <- rbind(combined,subset)}

      if (exists("factorgroupings.genes") == FALSE) {
        if(gaps.row==TRUE){
          ind.row <- ind.row + nrow(subset)
          gaps.row.n <- c(gaps.row.n,ind.row)
        }
      }

      }

    if (is.null(order.by.sample)==TRUE){combined <- combined[gene.order,]} #this is setting it, the order is determined in the loop, kind of redundant but dont feel like changing now

    if (sum(groupings == FALSE, na.rm=T) !=0) {  ##if not grouped by samples, see if should be ordered by gene, set clustering of samples
      if (is.null(order.by.gene)==FALSE){
        if(is.raw.Ct==FALSE){combined <- combined[,order(combined[which(rownames(combined) %in% order.by.gene),],na.last = F)]}
        if(is.raw.Ct==TRUE){combined <- combined[,order(combined[which(rownames(combined) %in% order.by.gene),],na.last = T)]}
        clust.cols <- F
      }else {
        #clustering of samples
        if (method %in% c("spearman","pearson", "kendall")) {
          clust.samps<-(as.dist(1-cor(combined,method=method,use=NA.handling)))
        }

        if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
          clust.samps <- dist(t(combined), method=method)
        }
      }

      gaps.col <- NULL
    }


    if(gaps.row==TRUE){
      if (exists("factorgroupings.genes") == FALSE) {gaps.row <- gaps.row.n[-length(gaps.row.n)]} else{gaps.row <- gaps.groupings.genes}
      gaps.row <- sort(rep(gaps.row,gap.width))}

    clust.rows <- F

    data.subset <- combined
  }




  #dealing with gaps and specific gaps of columns -- need to be able to set this for all scenarios

  if(is.null(gaps.col.spec)==FALSE){gaps.col <- gaps.col.spec; gaps.col <- sort(rep(gaps.col, gap.width))}
  if(is.null(gaps.row.spec)==FALSE){gaps.row <- gaps.row.spec; gaps.row <- sort(rep(gaps.row, gap.width))}



  if(clust.rows==T){heightrow <- treeheight.row}
  if(clust.cols==T){heightcol <- treeheight.col}

  #fixing colors and breaks at the end, regardless of separations --- not going to be combined anymore
  data.subset1 <- data.subset
  data.subset <- scales::squish(data.subset,params$scale.range)
  breaks <- seq(params$scale.range[1], params$scale.range[2],length.out=params$n.colors.range)
  my_cols=colorRampPalette(params$scale.colors)(n=params$n.colors.range-1)
  if(is.raw.Ct==TRUE){my_cols <- rev(my_cols)}

  #na.fix, regardless of separations
  if(na.fix==TRUE){
    if(is.raw.Ct==TRUE){
      data.subset[which(data.subset1==max(data.subset1))] <- params$scale.range[2]+0.04
      breaks <- c(breaks, params$scale.range[2]+0.01, params$scale.range[2]+0.05)
      my_cols <- c(my_cols,params$scale.colors[1],"grey90")
    }
    if(is.raw.Ct==FALSE){
      data.subset[which(data.subset1==min(data.subset1))] <- params$scale.range[1]-0.04
      breaks <- c(params$scale.range[1]-0.05,params$scale.range[1]-0.01,breaks)
      my_cols <- c("grey90",params$scale.colors[1],my_cols)
    }
  }

  temp.annot_samps <- params$annot_samps
  temp.annot_genes <- params$annot_genes
  temp.annot_cols <- params$annot_cols

  if (drop.annot.levels == TRUE) {
    suppressWarnings( if (is.na(temp.annot_samps) == F) {
      temp.annot_samps[] <- lapply(temp.annot_samps, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_samps <- temp.annot_samps %>% tibble::rownames_to_column("Sample")
      temp.annot_samps <- droplevels(temp.annot_samps[which(temp.annot_samps$Sample %in% colnames(data.subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="Sample")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_samps) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0 ) {
        for (annot.i in 1:length(spec.cols)) {
          annot <- colnames(temp.annot_samps)[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_samps[,which(colnames(temp.annot_samps)==annot)])  )]
          if ( sum( levels(temp.annot_samps[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_samps[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_samps[,annot])[levels(temp.annot_samps[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    })


    suppressWarnings( if (is.na(params$annot_genes) == F) {
      temp.annot_genes[] <- lapply(temp.annot_genes, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_genes <- temp.annot_genes %>% tibble::rownames_to_column("Gene")
      temp.annot_genes <- droplevels(temp.annot_genes[which(temp.annot_genes$Gene %in% rownames(data.subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="Gene")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_genes) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0) {
        for (annot.i in 1:length(colnames(temp.annot_genes))) {
          annot <- colnames(temp.annot_genes)[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_genes[,which(colnames(temp.annot_genes)==annot)])  )]
          if ( sum( levels(temp.annot_genes[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_genes[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_genes[,annot])[levels(temp.annot_genes[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    })

  }



  if (clust.cols == T) {
    tryclustcols <- try(hclust(clust.samps, linkage), silent = T)
    if (class(tryclustcols) == "try-error") {stop('cannot cluster columns, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.cols = F or specify order.by.gene')}
    }

  if (clust.rows == T) {
    tryclustrows <- try(hclust(clust.genes, linkage), silent = T)
    if (class(tryclustrows) == "try-error") {stop('cannot cluster rows, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.rows = F or specify order.by.sample')}
    }


  pheatmap(data.subset,col=my_cols, breaks=breaks, border_color = NA, clustering_method=linkage,annotation_col=temp.annot_samps, annotation_colors = temp.annot_cols,
           clustering_distance_rows = clust.genes, clustering_distance_cols = clust.samps, main=paste(main,"\n Method_",method,"_Linkage_",linkage),
           cluster_rows = clust.rows, cluster_cols = clust.cols, cutree_rows = row.groups, cutree_cols = col.groups, gaps_row = gaps.row, gaps_col = gaps.col,fontsize_row = fontsize.row, fontsize_col = fontsize.col,
           cellwidth = cell.width, cellheight = cell.height, show_rownames = show.rownames,show_colnames = show.colnames,
           treeheight_row = heightrow ,treeheight_col = heightcol, silent=hide.plot, legend=show.legend, annotation_legend = show.annotations,
           annotation_row=temp.annot_genes, drop_levels = drop.annot.levels)

    }


ExtractMatrix <- function(   ##from clustered heatmap, will extract the exact matrix
  data, ##data heatmap is based off of, (needs to be subset for genes already?)
  heatmap, ##input should be the heatmap saved as a variable, output will be the matrix in the order that the heatmap is output as
  clustered.cols = TRUE, ##assumes that heatmap$tree_col exists
  clustered.rows = TRUE ##assumes that heatmap$tree_row exists
){

  attributes <- unlist(lapply(heatmap$gtable$grobs, function(x)(x$name)))
  text.attributes <- grep("GRID.text", attributes)

  ordered <- data
  if(clustered.cols==TRUE){
    colord <- heatmap$gtable[["grobs"]][[text.attributes[2]]][["label"]]
    ordered <- ordered[,match(colord,colnames(ordered))]
  }

  if(clustered.rows==TRUE){
    roword <- heatmap$gtable[["grobs"]][[text.attributes[3]]][["label"]]
    ordered <- ordered[match(roword,rownames(ordered)),]
  }

  return(ordered)
}

extractClusters <- function(  ##will extract cluster groups from clustered heatmaps
  data, ##should be the same data as what was used to generate the heatmap, used to extract the correct gene order
  heatmap, #this should be the output from pheatmap that was saved
  extractGenes = FALSE, ##will extract from the gene clustering
  GeneGroup_Name = NULL,
  extractSamples = TRUE, ##will extract from the sample clustering
  SampleGroup_Name = NULL,
  num_Gene.groups = 5, ##default of 5, should be decided
  num_Sample.groups = 5, ##default of 5, should be decided
  all.genes = NULL, #gene names, if annotations should be in placed in larger gene list than matrix in question
  all.samples = NULL #sample names, if annotations should be in placed in larger sample list than matrix in question
){
  if(is.null(all.genes)==TRUE){ all.genes <- rownames(data) }
  if(is.null(all.samples)==TRUE){ all.samples <- colnames(data) }

  if (extractGenes == TRUE) {
    report.gene.groups <- rep("_",length(all.genes))
    gene.groups <- cutree(heatmap$tree_row, k=num_Gene.groups)
    alph <- LETTERS[1:num_Gene.groups]
    gene.tree.order <- heatmap$tree_row[["labels"]][heatmap$tree_row[["order"]]]
    ind=1
    for (i in 1:length(unique(gene.groups))) {
      which.group <- gene.groups[which(names(gene.groups)==gene.tree.order[ind])]
      names <- names(which(gene.groups==which.group))
      report.gene.groups[which(all.genes %in% names)] <- alph[i]
      ind <- ind + length(names)
    }
    report.gene.groups <- data.frame(Gene.Groups=report.gene.groups, stringsAsFactors = T); rownames(report.gene.groups) <- all.genes
    if (is.null(GeneGroup_Name) == FALSE) { colnames(report.gene.groups) <- GeneGroup_Name}
    #names(report.gene.groups) <- all.genes
  }

  if(extractSamples == TRUE){
    report.sample.groups <- rep("_",length(all.samples))
    sample.groups <- cutree(heatmap$tree_col,k=num_Sample.groups)
    alph <- LETTERS[1:num_Sample.groups]
    sample.tree.order <- heatmap$tree_col[["labels"]][heatmap$tree_col[["order"]]]
    ind=1
    for (i in 1:length(unique(sample.groups))) {
      which.group <- sample.groups[which(names(sample.groups)==sample.tree.order[ind])]
      names <- names(which(sample.groups==which.group))
      report.sample.groups[which(all.samples %in% names)] <- alph[i]
      ind <- ind + length(names)
    }
    report.sample.groups <- data.frame(Sample.Groups=report.sample.groups, stringsAsFactors = T); rownames(report.sample.groups) <- all.samples
    if (is.null(SampleGroup_Name) == FALSE) { colnames(report.sample.groups) <- SampleGroup_Name}
    # names(report.sample.groups) <- all.samples
  }

  if(extractGenes == TRUE & extractSamples == TRUE){return(c(list(Gene_Groups=report.gene.groups, Sample_Groups=report.sample.groups)))}
  if(extractGenes == TRUE & extractSamples == FALSE){return(report.gene.groups)}
  if(extractGenes == FALSE & extractSamples == TRUE){return(report.sample.groups)}
  #c(list(Gene_Groups=report.gene.groups, Sample_Groups=report.sample.groups))

}



extractGaps <- function(  ##similar to extractClusters, will return the positions of gaps created by separating cluster groups
  data, ##should be the same data as what was used to generate the heatmap, used to extract the correct gene order
  heatmap, #this should be the output from pheatmap that was saved
  extractRows = FALSE, ##will extract from the gene clustering
  extractCols = TRUE, ##will extract from the sample clustering
  num_Rows = 5, ##default of 5, should be decided
  num_Cols = 5 ##default of 5, should be decided
){
  if (extractRows == TRUE) {
    gene.groups <- cutree(heatmap$tree_row, k=num_Rows)
    gene.tree.order <- heatmap$tree_row[["labels"]][heatmap$tree_row[["order"]]]
    rowvec <- c(0)
    ind=1
    for (i in 1:length(unique(gene.groups))) {
      which.group <- gene.groups[which(names(gene.groups)==gene.tree.order[ind])]
      names <- names(which(gene.groups==which.group))
      rowvec <- c(rowvec,rowvec[i]+length(names))
      ind <- ind + length(names)
    }
    rowvec <- rowvec[-length(rowvec)]; rowvec <- rowvec[-1]
  }

  if(extractCols == TRUE){
    sample.groups <- cutree(heatmap$tree_col,k=num_Cols)
    sample.tree.order <- heatmap$tree_col[["labels"]][heatmap$tree_col[["order"]]]
    colvec <- c(0)
    ind=1
    for (i in 1:length(unique(sample.groups))) {
      which.group <- sample.groups[which(names(sample.groups)==sample.tree.order[ind])]
      names <- names(which(sample.groups==which.group))
      colvec <- c(colvec,colvec[i]+length(names))
      ind <- ind + length(names)
    }
    colvec <- colvec[-length(colvec)]; colvec <- colvec[-1]
  }

  if(extractRows == TRUE & extractCols == TRUE){return(c(list(Row_Groups=rowvec, Col_Groups=colvec)))}
  if(extractRows == TRUE & extractCols == FALSE){return(rowvec)}
  if(extractRows == FALSE & extractCols == TRUE){return(colvec)}
}



AOV1way <- function(
  data.to.aov,
  category,
  pthreshold = 0.05,
  additional.report = "NONE"  ##options are "NONE", "TUKEY","AOV", or "ALL"
){

  ####if data is not all samples, subset annotations appropriately
    if (sum(colnames(data.to.aov) %notin% rownames(params$annotations)) != 0 ) {
      stop('colnames of input data do not match rownames of annotations, cannot link annotations to data and assign groupings for ANOVA')}

  temp.annotations <- params$annotations[match(colnames(data.to.aov), rownames(params$annotations)),]

  groupings <- as.factor(droplevels(temp.annotations[,category]))

  aov.all <- apply(data.to.aov, 1, function(x)(summary(aov(x~groupings))))

  aov.results <- data.frame(FVal=unlist(lapply(aov.all,function(x)((x[[1]]$`F value`[1])))),
                            pVal=unlist(lapply(aov.all,function(x)((x[[1]]$`Pr(>F)`[1])))),
                            row.names = names(aov.all))

  sig.genes <- rownames(aov.results)[which(aov.results$pVal <= pthreshold)]
  nonsig.genes <- rownames(aov.results)[which(aov.results$pVal > pthreshold)]

  sig.set <- data.to.aov[which(rownames(data.to.aov) %in% sig.genes),]

  tukey.all <-  apply(sig.set,1,function(x)(TukeyHSD(aov(x~groupings))))

  tukey.pvals <- data.frame(Reduce(rbind, lapply(tukey.all,function(x)(x$groupings[,4]))), row.names = names(tukey.all)); colnames(tukey.pvals) <- gsub("\\.","-", colnames(tukey.pvals))
  tukey.diffs <- data.frame(Reduce(rbind, lapply(tukey.all,function(x)(x$groupings[,1]))), row.names = names(tukey.all)); colnames(tukey.diffs) <- gsub("\\.","-", colnames(tukey.diffs))



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



AOV2way <- function(
  data.to.aov,
  category1,
  category2,
  pthreshold = 0.05,
  additional.report = "NONE"  ##options are "NONE", "TUKEY","AOV", or "ALL"
){

  ####if data is not all samples, subset annotations appropriately
  if (sum(colnames(data.to.aov) %notin% rownames(params$annotations)) != 0 ) {
    stop('colnames of input data do not match rownames of annotations, cannot link annotations to data and assign groupings for ANOVA')}

  temp.annotations <- params$annotations[match(colnames(data.to.aov), rownames(params$annotations)),]

  groupings1 <- as.factor(droplevels(temp.annotations[,category1]))
  groupings2 <- as.factor(droplevels(temp.annotations[,category2]))

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

  tukey.pvals1 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% category1.sig)],function(x)(x$groupings1[,4]))), row.names = names(tukey.all)[which(names(tukey.all) %in% category1.sig)]); colnames(tukey.pvals1) <- rownames(tukey.all[[1]]$groupings1) # gsub("\\.","-", colnames(tukey.pvals1))
  tukey.pvals2 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% category2.sig)],function(x)(x$groupings2[,4]))), row.names = names(tukey.all)[which(names(tukey.all) %in% category2.sig)]); colnames(tukey.pvals2) <- rownames(tukey.all[[1]]$groupings2) #  gsub("\\.","-", colnames(tukey.pvals2))
  tukey.pvals3 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% interaction.sig)],function(x)(x$`groupings1:groupings2`[,4]))), row.names = names(tukey.all)[which(names(tukey.all) %in% interaction.sig)]); colnames(tukey.pvals3) <-  rownames(tukey.all[[1]]$`groupings1:groupings2`) # gsub("\\.","-", colnames(tukey.pvals3))

  tukey.diffs1 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% category1.sig)],function(x)(x$groupings1[,1]))), row.names = names(tukey.all)[which(names(tukey.all) %in% category1.sig)]); colnames(tukey.diffs1) <- rownames(tukey.all[[1]]$groupings1) # gsub("\\.","-", colnames(tukey.diffs1))
  tukey.diffs2 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% category2.sig)],function(x)(x$groupings2[,1]))), row.names = names(tukey.all)[which(names(tukey.all) %in% category2.sig)]); colnames(tukey.diffs2) <- rownames(tukey.all[[1]]$groupings2) #  gsub("\\.","-", colnames(tukey.diffs2))
  tukey.diffs3 <- data.frame(Reduce(rbind, lapply(tukey.all[which(names(tukey.all) %in% interaction.sig)],function(x)(x$`groupings1:groupings2`[,1]))), row.names = names(tukey.all)[which(names(tukey.all) %in% interaction.sig)]); colnames(tukey.diffs3) <-  rownames(tukey.all[[1]]$`groupings1:groupings2`) # gsub("\\.","-", colnames(tukey.diffs3))

  if (toupper(additional.report) == "ALL") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                "Category1-Sig.Genes" = category1.sig,
                "Category2-Sig.Genes" = category2.sig,
                "Interaction-Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Tukey.output' = tukey.all,
                'Category1-Tukey.pVals' = tukey.pvals1,
                'Category2-Tukey.pVals' = tukey.pvals2,
                'Interaction-Tukey.pVals' = tukey.pvals3,
                'Category1-Tukey.diffs' = tukey.diffs1,
                'Category2-Tukey.diffs' = tukey.diffs2,
                'Interaction-Tukey.diffs' = tukey.diffs3
    ))
  }


  if (toupper(additional.report) == "AOV") {
    return(list('AOV.output' = aov.all,
                'AOV.Results' = aov.results,
                "Category1-Sig.Genes" = category1.sig,
                "Category2-Sig.Genes" = category2.sig,
                "Interaction-Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1-Tukey.pVals' = tukey.pvals1,
                'Category2-Tukey.pVals' = tukey.pvals2,
                'Interaction-Tukey.pVals' = tukey.pvals3,
                'Category1-Tukey.diffs' = tukey.diffs1,
                'Category2-Tukey.diffs' = tukey.diffs2,
                'Interaction-Tukey.diffs' = tukey.diffs3
    ))
  }

  if (toupper(additional.report) == "TUKEY") {
    return(list('AOV.Results' = aov.results,
                "Category1-Sig.Genes" = category1.sig,
                "Category2-Sig.Genes" = category2.sig,
                "Interaction-Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1-Tukey.pVals' = tukey.pvals1,
                'Category2-Tukey.pVals' = tukey.pvals2,
                'Interaction-Tukey.pVals' = tukey.pvals3,
                'Category1-Tukey.diffs' = tukey.diffs1,
                'Category2-Tukey.diffs' = tukey.diffs2,
                'Interaction-Tukey.diffs' = tukey.diffs3
    ))
  }

  if (toupper(additional.report) == "NONE") {
    return(list('AOV.Results' = aov.results,
                "Category1-Sig.Genes" = category1.sig,
                "Category2-Sig.Genes" = category2.sig,
                "Interaction-Sig.Genes" = interaction.sig,
                "All.Sig.Genes" = any.sig,
                'NonSig.Genes' = nonsig.genes,
                'Category1-Tukey.pVals' = tukey.pvals1,
                'Category2-Tukey.pVals' = tukey.pvals2,
                'Interaction-Tukey.pVals' = tukey.pvals3,
                'Category1-Tukey.diffs' = tukey.diffs1,
                'Category2-Tukey.diffs' = tukey.diffs2,
                'Interaction-Tukey.diffs' = tukey.diffs3
    ))
  }
}




myPCA <- function(
  data,
  to.pca = "samples",
  nPcs= 3,
  color.by = "blue", ##vetor same length or in annotations, annot_cols
  custom.color.vec = FALSE,
  PCs.to.plot = c("PC1","PC2"),
  legend.position = "right",
  main = NULL,
  percent.mad =0.5
){

  if (to.pca == "samples") {

    temp.annotations <- params$annotations
    temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]


    pca <- pcaMethods::pca(t(data), nPcs = nPcs)
    pca.scrs <- pcaMethods::scores(pca)
    pca.ldgs <- pcaMethods::loadings(pca)

    pca.data <- data.frame(pca.scrs,temp.annotations)

    if (color.by %in% rownames(data) | sum(custom.color.vec != FALSE) > 0) {
      if (color.by %in% rownames(data)) {
        genedat<- data[which(rownames(data)==color.by),]
        cols <- myColorRamp5(params$expression_gradient.colors,genedat, percent.mad = percent.mad)
      } else{ cols <- custom.color.vec}

      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols))+ geom_point(pch=21,color="black",size=5)  +
        scale_fill_identity() +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)

    }else{


      if (color.by %in% colnames(temp.annotations)) {
        suppressWarnings( if (is.na(temp.annotations) == FALSE) {
          if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
            stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')}

          temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),]
        })

        if (color.by %in% names(params$annot_cols)) {
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
        }else{
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- hue_pal()(length(levels(cols)))
        }

      } else{ cols <- color.by; colors <- color.by}


      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols))+ geom_point(pch=21,color="black",size=5)  +
        scale_fill_manual(values=colors) +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
    }
  }

  if (to.pca == "genes") {

    temp.annotations.genes <- params$annotations.genes
    temp.annotations.genes <- temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),, drop = FALSE]


    pca <- pcaMethods::pca((data), nPcs = nPcs)
    pca.scrs <- pcaMethods::scores(pca)
    pca.ldgs <- pcaMethods::loadings(pca)

    pca.data <- data.frame(pca.scrs,temp.annotations.genes)

    if (sum(custom.color.vec != FALSE) > 0) {
      cols <- custom.color.vec

      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols))+ geom_point(pch=21,color="black",size=5)  +
        scale_fill_identity() +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)

    }else{


      if (color.by %in% colnames(temp.annotations.genes)) {
        suppressWarnings( if (is.na(temp.annotations.genes) == FALSE) {
          if (sum(rownames(data) %notin% rownames(temp.annotations.genes)) != 0 ) {
            stop('rownames of input data do not match rownames of annotations, cannot link annotations to data')}

          temp.annotations.genes <-temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),]
        })
        if (color.by %in% names(params$annot_cols)) {
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
        }else{
          cols <- as.factor(pca.data[,which(colnames(pca.data) == color.by)])
          colors <- hue_pal()(length(levels(cols)))
        }

      } else{ cols <- color.by; colors <- color.by}


      p <- ggplot(pca.data, aes(x=eval(parse(text = PCs.to.plot[1])),y=eval(parse(text = PCs.to.plot[2])),fill=cols))+ geom_point(pch=21,color="black",size=5)  +
        scale_fill_manual(values=colors) +labs(x=paste(PCs.to.plot[1]), y= paste(PCs.to.plot[2])) + ggtitle(main) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
    }
  }

  return(p)
}



PTM <- function(  #pavlidis template matching, correlation to a chosen template
  data, ##dataset, genes in the rows
  match.template, ##from params$annotations or gene, can also be gene or sample
  set.high, ##level(s) of annotations vector to set high
  custom.template = NA,  ##must be same length as data
  Find.Match.For = "genes", ##default to compare to a gene or sample metadata template, can change to "samples"
  cutoff = 0.05, ##default of 0.05, can be changed, assumes pval
  cut.by = "pvals", ##pvals will return lower than cutoff with a positive correlation, can also use rvals, returns values greater than cutoff
  method = "pearson", ##default of pearson
  NA.handling = "pairwise.complete.obs", ##default of pairwise complete obs, can be changed
  return.vals = FALSE
){

  if (tolower(Find.Match.For) == "genes") {   ##length of samples

    temp.annotations <- params$annotations



    if(is.na(custom.template)==TRUE){
      if (match.template %in% colnames(temp.annotations)) {
        suppressWarnings( if (is.na(temp.annotations) == FALSE) {
          if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
            stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
          }
          temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]}  )

        template <- rep(0,nrow(temp.annotations))
        template[which(temp.annotations[,match.template] %in% set.high)] <- 1
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

    suppressWarnings( if (is.na(temp.annotations.genes) == FALSE) {
      if (sum(rownames(data) %notin% rownames(temp.annotations.genes)) != 0 ) {
        stop('rownames of input data do not match rownames of annotations, cannot link annotations to data')
      }
      temp.annotations.genes <- temp.annotations.genes[match(rownames(data), rownames(temp.annotations.genes)),, drop = FALSE]}  )


    if(is.na(custom.template)==TRUE){
      if (match.template %in% colnames(temp.annotations.genes)) {
        template <- rep(0,nrow(temp.annotations.genes))
        template[which(temp.annotations.genes[,match.template] %in% set.high)] <- 1
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
              show.annotations = show.annotations, is.raw.Ct = is.raw.Ct, drop.annot.levels = drop.annot.levels)

    if(show.report==T){ report <- gen.cor[,which(gen.cor>above| gen.cor< below)]; return(report)}

  }

}


correlateGenes <- function(  ##broad gene correlations
  data, ##data matrix with genes in the rows
  limits =NULL, ## vector of two numbers, giving the lower and upper bounds
  nbreaks=20,
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  NA.handling = "pairwise.complete.obs"   ##use for correlations, can be overwritten
){
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



reportGenes <- function(  ##returns report summary of range of expression
  data, ##dataset, genes should be in rows
  list, ##list of genes to summarize
  ranges="fixed", ##default, can be changed, other option is mad
  fixed.range =2, #only if ranges=fixed, can be changed to adjust range
  weight = 1.25 #only used if ranges=mad, weight to be added to mad for ranges
){
  report <- matrix(nrow=length(list),ncol=6)

  if (sum(list %in% rownames(data)) == 0 ) {stop('matches for list not found in rownames data')}

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

  colnames(report) <- c("Gene","Percent Samples Detected","Percent Samples NA","Percent High Expressing",
                        "Percent Mid-Range Expression","Percent Low Expressing")

  return(as.data.frame(report))
}




makefactorgroup <- function(
  annots,
  levels,
  specify.gaps = NULL,
  return.gaps = FALSE
){

  if (length(levels) == 1) {
    Level.1 <- annots[,which(colnames(annots)==levels[1])]
    factor.group <- factor(Level.1); names(factor.group) <- rownames(annots)


    if (is.null(specify.gaps) == FALSE) {
      if (length(specify.gaps)  != length(levels)) {stop('length of gap specifications not equal to number of levels provided')}

      gaps <- rep(cumsum(rev(rev(rle(as.vector(sort(factor.group)))$lengths)[-1])), specify.gaps[1])

    }else{gaps <- cumsum(rev(rev(rle(as.vector(sort(factor.group)))$lengths)[-1]))}

  }

  if (length(levels) == 2) {
    Level.1 <- annots[,which(colnames(annots)==levels[1])]
    Level.2 <- annots[,which(colnames(annots)==levels[2])]
    combo <- data.frame((Level.1), (Level.2), paste(Level.1, Level.2)); colnames(combo) <- c("Level1","Level2","Combo"); rownames(combo) <- rownames(annots)
    bylevel1 <- combo[order(combo[,1]),]
    bylevel2 <- bylevel1[order(bylevel1[,2]),]

    factor.group <- factor(combo[,which(colnames(combo)=="Combo")], levels=c(unique(as.character(bylevel2[,which(colnames(bylevel2)=="Combo")])))); names(factor.group) <- rownames(annots)

    if (is.null(specify.gaps) == FALSE) {
      if (length(specify.gaps)  != length(levels)) {stop('length of gap specifications not equal to number of levels provided')}

      gaps <- sort(c( rep(cumsum(rev(rev(rle(as.vector(bylevel2$Level1))$lengths)[-1])), specify.gaps[1]),
                      rep(cumsum(rev(rev(rle(as.vector(bylevel2$Level2))$lengths)[-1])), specify.gaps[2])))

    }else{gaps <- cumsum(rev(rev(rle(as.vector(bylevel2$Combo))$lengths)[-1]))}
  }

  if (length(levels) == 3 ) {
    Level.1 <- annots[,which(colnames(annots)==levels[1])]
    Level.2 <- annots[,which(colnames(annots)==levels[2])]
    Level.3 <- annots[,which(colnames(annots)==levels[3])]

    combo <- data.frame((Level.1), (Level.2),(Level.3), paste(Level.1, Level.2,Level.3)); colnames(combo) <- c("Level1","Level2","Level3","Combo"); rownames(combo) <- rownames(annots)
    bylevel1 <- combo[order(combo[,1]),]
    bylevel2 <- bylevel1[order(bylevel1[,2]),]
    bylevel3 <- bylevel2[order(bylevel2[,3]),]

    factor.group <- factor(combo[,which(colnames(combo)=="Combo")], levels=c(unique(as.character(bylevel3[,which(colnames(bylevel3)=="Combo")])))); names(factor.group) <- rownames(annots)

    if (is.null(specify.gaps) == FALSE) {
      if (length(specify.gaps)  != length(levels)) {stop('length of gap specifications not equal to number of levels provided')}

      gaps <- sort(c( rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level1))$lengths)[-1])), specify.gaps[1]),
                      rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level2))$lengths)[-1])), specify.gaps[2]),
                      rep(cumsum(rev(rev(rle(as.vector(bylevel3$Level3))$lengths)[-1])), specify.gaps[3]) ) )

    }else{gaps <- cumsum(rev(rev(rle(as.vector(bylevel3$Combo))$lengths)[-1]))}

  }

  if (return.gaps == FALSE) {return(factor.group)}

  if (return.gaps == TRUE) {return(list(factor.group = factor.group, gaps = gaps))}

}




find.silhouette <- function(
  data,
  ngroups, #for rand to clust
  maxgroups=12,
  max.iter=10,
  method = "pearson",
  NA.handling = "pairwise.complete.obs",
  linkage = "complete",
  to.sil = "samples", # can change to genes
  to.view="rand.to.clust", ##can also be "all.clusts" or "rand.all.clusts"
  main="Average Silhouette Width",
  axis.label="Silhouette Width", #x axis for rand.to clust, y label for the others
  main.label.size= 30,
  axis.label.size=20,
  legend.position = "bottom"
){


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

    sil.coef <- summary(silhouette(groupings, clust.mat))$avg.width

    rands <-  replicate(max.iter,summary(silhouette(permute(groupings), clust.mat))$avg.width, simplify = "vector")

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
      sil.coefs <- c(sil.coefs, summary(silhouette(groupings, clust.mat))$avg.width)

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
        rands <- cbind(rands,replicate(max.iter,summary(silhouette(permute(groupings), clust.mat))$avg.width, simplify = "vector"))

      }


      to.plot <- data.frame(clusts=2:maxgroups, coefs=sil.coefs, type="Original Data")
      colnames(rands) <- 2:maxgroups
      rands.melt <- melt(rands); rands.melt$type <- "Randomized Data"



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


scatterGenes <- function(
  data,
  gene1,
  gene2,
  is.raw.Ct = FALSE, ##if data is raw and axis should be flipped, set to TRUE
  na.fix = 2,
  color.by= "blue",  ##can be a color, gene, otherwise utilize annot_samps and annot_cols
  custom.color.vec = FALSE, ##give custom vector, same order as samples
  xlimits = FALSE, #will make limits automatically, can switch to specify
  ylimits = FALSE,
  squish1 = FALSE, #if limits are specified, will remove points outside the range, can change to set to mins, maxs
  squish2 = FALSE,
  point.size = 5,
  transparency = 1,
  legend.position = "none",
  percent.mad = 0.5
){

  if (gene1 %notin% rownames(data)) {stop('gene1 not found in rownames data')}
  if (gene2 %notin% rownames(data)) {stop('gene2 not found in rownames data')}

  dat1<-data[which(rownames(data) %in% gene1),]; if (is.raw.Ct==F & na.fix!=F) {dat1[which(is.na(dat1))] <- (min(dat1, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat1[which(is.na(dat1))]<- (max(dat1, na.rm=T)+na.fix)}
  dat2<-data[which(rownames(data) %in% gene2),]; if (is.raw.Ct==F & na.fix!=F) {dat2[which(is.na(dat2))] <- (min(dat2, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat2[which(is.na(dat2))]<- (max(dat2, na.rm=T)+na.fix)}

  temp.annotations <- params$annotations
  temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]


  dat.to.plot <- data.frame(Gene1= dat1, Gene2= dat2); dat.to.plot <- cbind(dat.to.plot, temp.annotations)


  if (color.by %in% rownames(data) | sum(custom.color.vec != FALSE) > 0) {
    if (color.by %in% rownames(data)) {
      genedat<- data[which(rownames(data)==color.by),]
      cols <- myColorRamp5(params$expression_gradient.colors,genedat, percent.mad = percent.mad)
    } else{ cols <- custom.color.vec}

    if (((xlimits==FALSE) && (ylimits==FALSE)) == TRUE) {
      suppressWarnings( if (squish1 != FALSE) {dat.to.plot$Gene1 <-  scales::squish(dat.to.plot$Gene1,squish1)} )
      suppressWarnings( if (squish2 != FALSE) {dat.to.plot$Gene2 <-  scales::squish(dat.to.plot$Gene2,squish2)} )

      if (is.raw.Ct==T) {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=5, alpha = transparency)  +
          scale_fill_identity() +labs(x=paste(gene1), y= paste(gene2)) +ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)+
          scale_x_reverse() + scale_y_reverse()
      } else {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=5, alpha = transparency)  +
          scale_fill_identity() +labs(x=paste(gene1), y= paste(gene2)) + ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
      }
    }


    if ((xlimits || ylimits) == TRUE) {
      # if (squish.to.lims != FALSE) {dat.to.plot$Gene1 <-  scales::squish(dat.to.plot$Gene1,xlimits)
      # dat.to.plot$Gene2 <-  scales::squish(dat.to.plot$Gene2,ylimits)}

      if (is.raw.Ct==T) {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=5, alpha = transparency)  +
          scale_fill_identity() +labs(x=paste(gene1), y= paste(gene2)) +ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)+
          xlim(c(xlimits)) + ylim(c(ylimits)) #+ scale_x_reverse() + scale_y_reverse()

      } else {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=5, alpha = transparency)  +
          scale_fill_identity() +labs(x=paste(gene1), y= paste(gene2)) + ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position) + xlim(c(xlimits)) + ylim(c(ylimits))
      }
    }
  } else{

    if (color.by %in% colnames(temp.annotations)) {
        if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
          stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
        }
        temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]

      if (color.by %in% names(params$annot_cols)) {
        cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
        colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
      }else{
        cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
        colors <- hue_pal()(length(levels(cols)))
      }
    } else{ cols <- color.by; colors <- color.by}


    if (((xlimits==FALSE) && (ylimits==FALSE)) == TRUE) {
      suppressWarnings( if (squish1 != FALSE) {dat.to.plot$Gene1 <-  scales::squish(dat.to.plot$Gene1,squish1)} )
      suppressWarnings( if (squish2 != FALSE) {dat.to.plot$Gene2 <-  scales::squish(dat.to.plot$Gene2,squish2)} )

      if (is.raw.Ct==T) {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
          scale_fill_manual(values=colors) +labs(x=paste(gene1), y= paste(gene2)) +ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)+
          scale_x_reverse() + scale_y_reverse()
      } else {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
          scale_fill_manual(values=colors) +labs(x=paste(gene1), y= paste(gene2)) + ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)
      }
    }


    if ((xlimits || ylimits) == TRUE) {
      # if (squish.to.lims != FALSE) {dat.to.plot$Gene1 <-  scales::squish(dat.to.plot$Gene1,xlimits)
      # dat.to.plot$Gene2 <-  scales::squish(dat.to.plot$Gene2,ylimits)}

      if (is.raw.Ct==T) {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
          scale_fill_manual(values=colors) +labs(x=paste(gene1), y= paste(gene2)) +ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position)+
          xlim(c(xlimits)) + ylim(c(ylimits)) #+ scale_x_reverse() + scale_y_reverse()

      } else {
        p <- ggplot(dat.to.plot, aes(x=Gene1,y=Gene2,fill=cols))+ geom_point(pch=21,color="black",size=point.size, alpha = transparency)  +
          scale_fill_manual(values=colors) +labs(x=paste(gene1), y= paste(gene2)) + ggtitle(paste(gene2, "vs.",gene1)) +
          theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                             axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = legend.position) + xlim(c(xlimits)) + ylim(c(ylimits))
      }
    }
  }

  return(p)
}


beeswarmGenes <- function( ##can save as ggplot object and add layers afterwards if more specifications need to be changed
  data,
  list,
  exact = TRUE,
  is.raw.Ct = FALSE,
  na.fix = 2,
  squishy = FALSE,      ##might need to add option for limits as well
  color.by = "blue",   ##single color, gene, column in annot_samps and use annot_cols
  custom.color.vec = FALSE,
  groupby.x = NULL, #option to change what is grouped by or on the X axis if faceted, change to false if groups arent needed,  #if null and color.by is not in annot_samps, will not facet and will not split into groups, equivalent to setting equal to FALSE
  custom.group.vec = FALSE,
  facet.wrap = FALSE, ##can change to true
  ncols=2, ##can change
  scales="free_y",
  legend.position = "none",
  axis.text.x.size = 25,
  point.size = 3,
  transparency = 1,
  percent.mad = 0.5
){

  ###set up, get genes, squish scale if needed, set groupby.x == FALSE if it doesnt match with colors
  if (exact == TRUE) {dat<-data[which(rownames(data) %in% list),]
    if (length(dat) == 0 ) {stop('exact matches for list not found in rownames data')}
    if (is.raw.Ct==F & na.fix!=F) {dat[which(is.na(dat))] <- (min(dat, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat[which(is.na(dat))]<- (max(dat, na.rm=T)+na.fix)}}
  if (exact == FALSE) {dat<-data[grep(paste(list, collapse = "|"),rownames(data)),]
    if (length(dat) == 0 ) {stop('inexact matches for list not found in rownames data')}
    if (is.raw.Ct==F & na.fix!=F) {dat[which(is.na(dat))] <- (min(dat, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat[which(is.na(dat))]<- (max(dat, na.rm=T)+na.fix)}}


  ####if data is not all samples, subset annotations appropriately
  temp.annotations <- params$annotations
  temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]

  suppressWarnings( if (is.na(temp.annotations) == FALSE) {
    if (sum(colnames(dat) %notin% rownames(temp.annotations)) != 0 ) {
      stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
    }
    temp.annotations <- temp.annotations[match(colnames(dat), rownames(temp.annotations)),, drop = FALSE]
  }  )

  if (is.null(groupby.x) == TRUE & (color.by %in% colnames(temp.annotations)) == FALSE) { groupby.x <- FALSE}  ##if groupby.x is null and color.by is in annot_samps, will group by that annotation as well, if no override to group and no annotation to color, wont group at all, if custom group vector supplied, will get corrected downstream


  ####if coloring by gene or custom color vector, identity based

  if (color.by %in% rownames(data) | sum(custom.color.vec != FALSE) > 0) {   ##if coloring by gene or by custom
    if (color.by %in% rownames(data)) {
      genedat<- data[which(rownames(data)==color.by),]
      if (is.raw.Ct ==FALSE) {cols <- myColorRamp5(params$expression_gradient.colors,genedat, percent.mad = percent.mad)}
      if (is.raw.Ct ==TRUE) {cols <- myColorRamp5(rev(params$expression_gradient.colors),genedat, percent.mad = percent.mad)}
    } else{ cols <- custom.color.vec}


    ##make dat.to.plot with identiy based colors
    suppressWarnings( if (custom.group.vec != FALSE) {
      dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations); dat.to.plot$cols <- cols; dat.to.plot$Custom <- custom.group.vec

      dat.to.plot <- melt(dat.to.plot, id.vars = c(colnames(temp.annotations),"cols", "Custom"))
      if (is.na(temp.annotations) == TRUE) {
        dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable == "temp.annotations"),]
      }

      groupby.x <- "Custom"

      suppressWarnings( if (squishy != FALSE) { dat.to.plot$value <- scales::squish(dat.to.plot$value, squishy)} )  ##if we want to squish
    }else{
      dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations); dat.to.plot$cols <- cols

      dat.to.plot <- melt(dat.to.plot, id.vars = c(colnames(temp.annotations),"cols"))
      if (is.na(temp.annotations) == TRUE) {
        dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable == "temp.annotations"),]
      }

      suppressWarnings( if (squishy != FALSE) { dat.to.plot$value <- scales::squish(dat.to.plot$value, squishy)} )  ##if we want to squish
    })   ##set dat.to.plot with identity based color vector and identity based group vector if supplied


    #####
    if ((is.null(groupby.x) == FALSE)) {  ##groupby has either been set to false by user or by previous tested condition (same as color, taken care of above)
      if (groupby.x != FALSE) {  ##set group to specification
        if (facet.wrap == FALSE) {
          if(is.raw.Ct==T){
            (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_identity() + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() )
          }else{
            (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_identity() + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") )
          }
        }

        if (facet.wrap == TRUE) {
          if(is.raw.Ct==T){
            (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=groupby.x)),y=value,fill=cols,group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, alpha = transparency) +
               scale_fill_identity() + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() +
               facet_wrap(~variable, ncol=ncols, scales = scales) )
          }else{
            (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=groupby.x)),y=value,fill=cols,group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size,alpha = transparency) +
               scale_fill_identity() + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") +
               facet_wrap(~variable, ncol=ncols, scales = scales) )

          }
        }
      }else{ ##no groupings, and no facet
        if(is.raw.Ct==T){
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", size=point.size) +
             scale_fill_identity() + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() )
        }else{
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", size=point.size, alpha = transparency) +
             scale_fill_identity() + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") )
        }
      }

    }

  }else{ ##custom color vector, by gene or supplied


    ####if color.by is by an annotation, not identity based colors

    ###set dat.to.plot
    suppressWarnings( if (custom.group.vec != FALSE) {
      dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations); dat.to.plot$Custom <- custom.group.vec

      dat.to.plot <- melt(dat.to.plot, id.vars = c(colnames(temp.annotations),"Custom"))
      if (is.na(temp.annotations) == TRUE) {
        dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable == "temp.annotations"),]
      }

      groupby.x <- "Custom"

      suppressWarnings( if (squishy != FALSE) { dat.to.plot$value <- scales::squish(dat.to.plot$value, squishy)} )  ##if we want to squish
    }else{
      dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations)

      dat.to.plot <- melt(dat.to.plot, id.vars = colnames(temp.annotations))
      if (is.na(temp.annotations) == TRUE) {
        dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable == "temp.annotations"),]
      }

      suppressWarnings( if (squishy != FALSE) { dat.to.plot$value <- scales::squish(dat.to.plot$value, squishy)} )  ##if we want to squish
    })

    if (color.by %in% colnames(temp.annotations)) {
      if (color.by %in% names(params$annot_cols)) {
        cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
        colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
      }else{
        cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
        colors <- hue_pal()(length(levels(cols)))
      }

    } else{ cols <- color.by; colors <- color.by} ##single color


    ##group by same annotations as coloring
    if ( (is.null(groupby.x) == TRUE) & (color.by %in% colnames(temp.annotations))) {

      if (facet.wrap == FALSE) {
        if(is.raw.Ct==T){
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=color.by))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() )
        }else{
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=color.by))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") )

        }
      }

      if (facet.wrap == TRUE) {
        if(is.raw.Ct==T){
          (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=color.by)),y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() +
             facet_wrap(~variable, ncol=ncols, scales = scales) )
        }else{
          (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=color.by)),y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") +
             facet_wrap(~variable, ncol=ncols, scales = scales) )

        }
      }

    }


    if ((is.null(groupby.x) == FALSE)) {
      if (groupby.x != FALSE) {  ##group by specified grouping
        if (facet.wrap == FALSE) {
          if(is.raw.Ct==T){
            (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() )
          }else{
            (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols, group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") )
          }
        }

        if (facet.wrap == TRUE) {
          if(is.raw.Ct==T){
            (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=groupby.x)),y=value,fill=cols,group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() +
               facet_wrap(~variable, ncol=ncols, scales = scales) )
          }else{
            (p <- ggplot(dat.to.plot, aes(x=eval(parse(text=groupby.x)),y=value,fill=cols,group=eval(parse(text=groupby.x))))+ geom_quasirandom(pch=21,color="black", dodge.width = 0.8, size=point.size, alpha = transparency) +
               scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
               theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                  strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                  axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") +
               facet_wrap(~variable, ncol=ncols, scales = scales) )

          }
        }
      }else{  ##set to false, no grouping and no faceting
        if(is.raw.Ct==T){
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Raw Ct Value") + scale_y_reverse() )
        }else{
          (p <- ggplot(dat.to.plot, aes(x=variable,y=value,fill=cols))+ geom_quasirandom(pch=21,color="black", size=point.size, alpha = transparency) +
             scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
             theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                                strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                                axis.title.y = element_text(size=20), axis.title.x=element_blank(), axis.text.x = element_text(size=axis.text.x.size)) + ylab("Normalized Expression Level") )
        }
      }

    }
  }
  return(p)
}



volcano <- function(
  data, ##dataset, genes should be in rows
  groups, ##vector the same length as number of samples, separating the two groups
  levels = NULL, ##levels of the groups, list first first, only list 2, if groups has more than two levels, pick levels
  is.log2 = TRUE, ##is data in log2 space? needed for FC vs LFC
  pval.cut =0.05, ##places horizontal line
  FC.cut= 2, ##five the fold change, function will put it in log2
  return.summary = FALSE,
  downreg.color = "green",
  upreg.color = "red",
  nosig.color = "gray",
  show.genes = NULL,
  point.size = 2,
  transparency = 1,
  legened.position = "right"
){

  ####if data is not all samples, subset annotations appropriately
  temp.annotations <- params$annotations


  if (groups %in% colnames(temp.annotations)) {

      if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
        stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
      }
      temp.annotations <- temp.annotations[match(colnames(data), rownames(temp.annotations)),, drop = FALSE]

    groupings <- as.factor(temp.annotations[,groups] )
    if (is.null(levels) == TRUE) { levels <- levels(groupings)}
    G1 <- data[,which(groupings==levels[1])]
    G2 <- data[,which(groupings==levels[2])]
  }else{
    G1 <- data[,which(groups==levels[1])]
    G2 <- data[,which(groups==levels[2])]
  }

  pvals <- NULL
  log2foldchanges <- NULL
  for (i in 1:nrow(G1)) {
    ttest <- t.test(G1[i,],G2[i,])
    pvals <- c(pvals,ttest$p.value)

    if (is.log2 == TRUE) {
      log2foldch <- ttest$estimate[2]-ttest$estimate[1]
    } else{ log2foldch <- log2(ttest$estimate[2]/ttest$estimate[1])}

    log2foldchanges <- c(log2foldchanges, log2foldch)

  }
  names(pvals) <- rownames(data); names(log2foldchanges) <- rownames(data)


  volcano.summary <- data.frame("LFC"=log2foldchanges,"FoldChange"=2^(log2foldchanges), pvals,"-log10pvals"=-log10(pvals))


  group <- rep("No Sig",nrow(volcano.summary))
  group[which(volcano.summary$pvals < pval.cut & (volcano.summary$LFC) > log2(FC.cut))] <-  "Upregulated" #paste("Fold Change", FC.cut, "& PValue <" pval.cut)  ##things that pass the original cutoff and p value
  group[which(volcano.summary$pvals < pval.cut & (volcano.summary$LFC) < -log2(FC.cut))] <- "Downregulated" #paste("Fold Change -", FC.cut, "& PValue <" pval.cut)  ##things that pass the original cutoff and p value

  mat <- cbind(volcano.summary,Color=group, Gene = rownames(volcano.summary))
  Sig.Genes <- rownames(volcano.summary); Sig.Genes[which(group == "No Sig")] <- ""
  mat <- cbind(mat, Sig.Genes)

  if (is.null(show.genes) == FALSE) {
    My.Genes <- rownames(volcano.summary); My.Genes[which(rownames(volcano.summary) %notin% show.genes)] <- ""
    mat <- cbind(mat, My.Genes)
  }

  p <- ggplot(mat,aes(x=LFC, y=-log10(pvals), col=Color)) + geom_point(size=2) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill="white"), panel.border = element_rect(color = "black", fill=NA), strip.background = element_blank(),
          strip.text = element_text(size=25), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title = element_text(size=20), plot.title = element_text(size=15, hjust = 0.5), legend.position = legend.position) +
    xlab("Log2 Fold Change") + ylab("-log10(Pvalue)") + scale_color_manual(name = paste(paste0("FC.cut = ", FC.cut), paste0("Pval.cut = ", pval.cut), sep="\n"), values=c("Downregulated"=downreg.color,"Upregulated"=upreg.color,"No Sig"=nosig.color)) +
    ggtitle(paste("-log10(pvalue) vs. log2(Fold Change) for",levels[2],"over",levels[1]))

  if (return.summary == FALSE) {return(p)}
  if (return.summary == TRUE) {return(volcano.summary)}

}


DensityGenes <- function(
  data,
  list,
  color.by = "blue", ##also dictates how it will split, need option to make custom vector to split on
  exact = TRUE,
  is.raw.Ct = FALSE,
  na.fix = 2,
  transparency = 0.5,
  ncols=2, ##can change
  scales="free",
  legend.position = "right"
){

  if (exact == TRUE) {dat<-data[which(rownames(data) %in% list),]
    if (length(dat) == 0 ) {stop('exact matches for list not found in rownames data')}
    if (is.raw.Ct==F & na.fix!=F) {dat[which(is.na(dat))] <- (min(dat, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat[which(is.na(dat))]<- (max(dat, na.rm=T)+na.fix)}}
  if (exact == FALSE) {dat<-data[grep(paste(list, collapse = "|"),rownames(data)),]
    if (length(dat) == 0 ) {stop('inexact matches for list not found in rownames data')}
    if (is.raw.Ct==F & na.fix!=F) {dat[which(is.na(dat))] <- (min(dat, na.rm=T)-na.fix)};if (is.raw.Ct==T & na.fix!=F) {dat[which(is.na(dat))]<- (max(dat, na.rm=T)+na.fix)}}


  temp.annotations <- params$annotations

  if (color.by %in% colnames(temp.annotations)) {

      if (sum(colnames(dat) %notin% rownames(temp.annotations)) != 0 ) {
        stop('colnames of input data do not match rownames of annotations, cannot link annotations to data')
      }
      temp.annotations <- temp.annotations[match(colnames(dat), rownames(temp.annotations)),, drop = FALSE]


    dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations)

    dat.to.plot <- melt(dat.to.plot, id.vars = colnames(temp.annotations))


    if (color.by %in% names(params$annot_cols)) {
      cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
      colors <- params$annot_cols[[which(names(params$annot_cols) == color.by)]]
    }else{
      cols <- as.factor(dat.to.plot[,which(colnames(dat.to.plot) == color.by)])
      colors <- hue_pal()(length(levels(cols)))
    }

    if(is.raw.Ct==T){
      p <- ggplot(dat.to.plot, aes(x=value,fill=cols, group=eval(parse(text = color.by))))+ geom_density(alpha = transparency) + facet_wrap(~variable, ncol=ncols, scales=scales) +
        scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                           axis.title = element_text(size=20)) + xlab("Raw Ct Value") + ylab("Denstiy") + scale_y_reverse()
    }else{
      p <- ggplot(dat.to.plot, aes(x=value,fill=cols, group=eval(parse(text = color.by))))+ geom_density(alpha = transparency) +  facet_wrap(~variable, ncol=ncols, scales=scales) +
        scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
        theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                           strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                           axis.title = element_text(size=20), axis.text.x = element_text(size=15)) + xlab("Normalized Expression Level") +ylab("Density")

    }
  } else{ cols <- color.by; colors <- color.by

  dat.to.plot <- data.frame(t(dat)); dat.to.plot <- cbind(dat.to.plot, temp.annotations)

  dat.to.plot <- melt(dat.to.plot, id.vars = colnames(temp.annotations))
  suppressWarnings(if (is.na(temp.annotations) == TRUE) {
    dat.to.plot <- dat.to.plot[-which(dat.to.plot$variable == "temp.annotations"),]
    })

  if(is.raw.Ct==T){
    p <- ggplot(dat.to.plot, aes(x=value,fill=cols))+ geom_density(alpha = transparency) + facet_wrap(~variable, ncol=ncols, scales=scales) +
      scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
      theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                         strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                         axis.title = element_text(size=20)) + xlab("Raw Ct Value") + ylab("Denstiy") + scale_y_reverse()
  }else{
    p <- ggplot(dat.to.plot, aes(x=value,fill=cols))+ geom_density(alpha = transparency) +  facet_wrap(~variable, ncol=ncols, scales=scales) +
      scale_fill_manual(values=colors) + #ggtitle(paste(list)) +
      theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40),
                         strip.text = element_text(size=25), strip.background.x = element_blank(), legend.position = legend.position,
                         axis.title = element_text(size=20), axis.text.x = element_text(size=15)) + xlab("Normalized Expression Level") +ylab("Density")

  }}
  return(p)
}

