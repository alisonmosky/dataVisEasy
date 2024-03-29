####Heatmap Functions####

#' @import pheatmap
#' @importFrom scales squish
#' @importFrom tibble column_to_rownames rownames_to_column remove_rownames
#' @importFrom stats sd
#' @export
myHeatmap <- function(  ##basic heatmap, can subset for gene list
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list = NULL,   ##list of gene or genes that will be extracted from the dataset, exact matches
  exact = TRUE, ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  linkage = "complete",   ##linkage method, defulat complete linkage, can be changed
  NA.handling = "pairwise.complete.obs",   ##use for correlations, can be overwritten
  clust.rows = TRUE, ##default for clustering rows, can be overwritten if error
  clust.cols = TRUE, ##same as clust.cols but for cols
  scale.rows = FALSE, ## option to zscore or median center rows
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
  show.rownames=TRUE,
  show.colnames=FALSE,
  treeheight.row=20,
  treeheight.col=20,
  hide.plot=FALSE,
  na.fix=FALSE,
  na.offset = 2,
  show.legend=TRUE,
  show.annotations=TRUE,
  is.raw.Ct=FALSE, ##if true, will reverse color scale to show yellow as high expressing
  drop.annot.levels=TRUE,
  annotation.names.row = TRUE,
  annotation.names.col = TRUE,
  border.color = NA,
  na.color = "grey90"
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }


  if (scale.rows != FALSE) {
    if (scale.rows %notin% c("median","zscore")) {stop('accetable methods for row scaling are "median" and "zscore"')
    }else{
      if (scale.rows == "zscore") {
        m <- apply(data, 1, mean, na.rm = T)
        s <- apply(data, 1, sd, na.rm = T)
        data <- ((data - m) / s)
      }
      if (scale.rows == "median") {
        data <- t(apply(data,1,function(x)(x-median(x, na.rm = T))))
      }
     }
  }


  if (is.null(list) == TRUE) {list <- rownames(data)}


  ##subset for list
  if (exact == TRUE) {
    subset <- as.matrix(data)[which(rownames(data) %in% list),]
    if (length(subset) == 0 ) {stop('exact matches for list not found in rownames data')}
  }else{
    subset <- as.matrix(data)[grep(paste(list, collapse = "|"),rownames(data)),]
    if (length(subset) == 0 ) {stop('inexact matches for list not found in rownames data')}
  }


  if (na.fix==TRUE) {
    if (is.raw.Ct==TRUE) {subset[which(is.na(subset))] <- max(subset,na.rm=T)+na.offset}
    if (is.raw.Ct==FALSE) {subset[which(is.na(subset))] <- min(subset,na.rm=T)-na.offset}
  }


  if (is.null(order.by.gene)==FALSE){
    if(is.raw.Ct==FALSE){subset <- subset[,order(data[which(rownames(data) %in% order.by.gene),],na.last = F)]}
    if(is.raw.Ct==TRUE){subset <- subset[,order(data[which(rownames(data) %in% order.by.gene),],na.last = T)]}
    clust.cols <- FALSE
  }

  if (is.null(order.by.sample)==FALSE){
    if(is.raw.Ct==FALSE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = F),]}
    if(is.raw.Ct==TRUE){subset <- subset[order(subset[,which(colnames(subset) %in% order.by.sample)],na.last = T),]}
    clust.rows <- FALSE
  }

  if(clust.rows==TRUE){heightrow <- treeheight.row}
  if(clust.cols==TRUE){heightcol <- treeheight.col}


  if(is.raw.Ct==TRUE){my_cols <- rev(my_cols)}

  if (any( rownames(params$annot_samps) %in% colnames(subset))) {
    temp.annot_samps <- params$annot_samps
  }else{ temp.annot_samps <- NA}

  if (any( rownames(params$annot_genes) %in% rownames(subset))) {
    temp.annot_genes <- params$annot_genes
  }else{ temp.annot_genes <- NA}

  temp.annot_cols <- params$annot_cols

  if (drop.annot.levels == TRUE) {
    if (any(!is.na(temp.annot_samps))) {
      temp.annot_samps[] <- lapply(temp.annot_samps, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_samps <- temp.annot_samps %>% tibble::rownames_to_column("rowN")
      temp.annot_samps <- droplevels(temp.annot_samps[which(temp.annot_samps$rowN %in% colnames(subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="rowN")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_samps) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0 ) {
        for (annot.i in 1:length(spec.cols)) {
          annot <- spec.cols[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_samps[,which(colnames(temp.annot_samps)==annot)])  )]
          if ( sum( levels(temp.annot_samps[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_samps[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_samps[,annot])[levels(temp.annot_samps[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    }


    if (any(!is.na(temp.annot_genes))) {
      temp.annot_genes[] <- lapply(temp.annot_genes, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_genes <- temp.annot_genes %>% tibble::rownames_to_column("rowN")
      temp.annot_genes <- droplevels(temp.annot_genes[which(temp.annot_genes$rowN %in% rownames(subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="rowN")

      spec.cols <- colnames(temp.annot_genes)[colnames(temp.annot_genes) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0) {
        for (annot.i in 1:length(spec.cols)) {
          annot <- spec.cols[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_genes[,which(colnames(temp.annot_genes)==annot)])  )]
          if ( sum( levels(temp.annot_genes[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_genes[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_genes[,annot])[levels(temp.annot_genes[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    }

  }




  if (clust.cols == TRUE) {
    if (method %in% c("spearman","pearson", "kendall")) {
      clust.samps<-(as.dist(1-cor(subset,method=method,use=NA.handling)))
    }
    if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
      clust.samps <- dist(t(subset), method=method)
    }

    tryclustcols <- try(hclust(clust.samps, linkage), silent = TRUE)
    if (class(tryclustcols) == "try-error") {stop('cannot cluster columns, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.cols = F or specify order.by.gene')}
  }else{
    if (is.null(gaps.col) == FALSE) {gaps.col <- sort(rep(gaps.col, gap.width))
    }
  }

  if (clust.rows == TRUE) {
    if (method %in% c("spearman","pearson", "kendall")) {
      clust.genes<-(as.dist(1-cor(t(subset),method=method,use=NA.handling)))
    }
    if (method %in% c("euclidean","maximum","manhattan","canberra","binary","minkowski")){
      clust.genes <- dist(subset,method=method)
    }

    tryclustrows <- try(hclust(clust.genes, linkage), silent = TRUE)
    if (class(tryclustrows) == "try-error") {stop('cannot cluster rows, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.rows = F or specify order.by.sample')}
  }else{
    if (is.null(gaps.row) == FALSE) {gaps.row <- sort(rep(gaps.row, gap.width))
    }
  }

  subset1 <- subset
  subset <- scales::squish(subset,params$scale.range)
  breaks <- seq(params$scale.range[1], params$scale.range[2],length.out=params$n.colors.range)
  my_cols=colorRampPalette(params$scale.colors)(n=params$n.colors.range-1)

  if(na.fix==TRUE){
    if(is.raw.Ct==TRUE){
      subset[which(subset1==max(subset1))] <- params$scale.range[2]+0.04
      breaks <- c(breaks, params$scale.range[2]+0.01, params$scale.range[2]+0.05)
      my_cols <- c(my_cols,params$scale.colors[1],na.color)      ##may need an option to set na_col
    }
    if(is.raw.Ct==FALSE){
      subset[which(subset1==min(subset1))] <- params$scale.range[1]-0.04
      breaks <- c(params$scale.range[1]-0.05,params$scale.range[1]-0.01,breaks)
      my_cols <- c(na.color,params$scale.colors[1],my_cols)
    }
  }

  if (is.null(main)==TRUE){
    main <- paste("Genes of Interest:",paste(list, collapse = ","))
    main <- paste(main,"\n Method:",method," Linkage:",linkage)}


  pheatmap(subset,col=my_cols, breaks=breaks, border_color = border.color, na_col = na.color, clustering_method=linkage,annotation_col=temp.annot_samps, annotation_colors = temp.annot_cols,
           clustering_distance_rows = clust.genes, clustering_distance_cols = clust.samps, main=main,
           cluster_rows = clust.rows, cluster_cols = clust.cols, cutree_rows = row.groups, cutree_cols = col.groups, gaps_row = gaps.row, gaps_col = gaps.col,
           cellwidth = cell.width, cellheight = cell.height, fontsize_row = fontsize.row, fontsize_col = fontsize.col, show_rownames = show.rownames,show_colnames = show.colnames,
           treeheight_row = heightrow ,treeheight_col = heightcol, silent = hide.plot, legend=show.legend, annotation_legend = show.annotations,
           annotation_row=temp.annot_genes, drop_levels = drop.annot.levels, annotation_names_row = annotation.names.row, annotation_names_col = annotation.names.col)


}




#' @export
myHeatmapByAnnotation <- function(
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list = NULL,   ##list of gene or genes that will be extracted from the dataset, exact matches
  exact = TRUE, ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
  groupings = FALSE, ##either character vector pointing to annotations, dataframe where the first row will be taken, or factor, if unnamed factor, wont properly order and will assume in same order as data
  groupings.gaps = NULL,
  groupings.genes = FALSE,
  groupings.genes.gaps = NULL,
  method = "pearson",   ##clustering and correlating method, default of pearson, can be switched to spearman
  linkage = "complete",   ##linkage method, defulat complete linkage, can be changed
  NA.handling = "pairwise.complete.obs",   ##use for correlations, can be overwritten
  clust.rows = TRUE, ##default for clustering rows, can be overwritten if error
  clust.cols = TRUE, ##same as clust.cols but for cols
  scale.rows = FALSE, ## option to zscore or median center rows
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
  drop.annot.levels = TRUE,
  annotation.names.row = TRUE,
  annotation.names.col = TRUE,
  border.color = NA,
  na.color = "grey90"
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }


  if (scale.rows != FALSE) {
    if (scale.rows %notin% c("median","zscore")) {stop('accetable methods for row scaling are "median" and "zscore"')
    }else{
      if (scale.rows == "zscore") {
        m <- apply(data, 1, mean, na.rm = T)
        s <- apply(data, 1, sd, na.rm = T)
        data <- ((data - m) / s)
      }
      if (scale.rows == "median") {
        data <- t(apply(data,1,function(x)(x-median(x, na.rm = T))))
      }
    }
  }


  if (is.null(list) == TRUE) {list <- rownames(data)}

  ##subset for list
  if (exact == TRUE) {
    data.subset <- as.matrix(data[which(rownames(data) %in% list),]);colnames(data.subset) <- colnames(data)
    if (length(data.subset) == 0 ) {stop('exact matches for list not found in rownames data')}

  }else{
    data.subset <- as.matrix(data[grep(paste(list, collapse = "|"),rownames(data)),]); colnames(data.subset) <- colnames(data)
    if (length(data.subset) == 0 ) {stop('exact matches for list not found in rownames data')}

  }


  ####new code to order annotations if not in order or if extra etc####
  if (any(!is.na(params$annotations))) {
    if (any(colnames(data.subset) %notin% rownames(params$annotations))) {
      warning('colnames of input data do not match rownames of annotations, cannot link annotations to data')
    }
    temp.annotations <- params$annotations[match(colnames(data.subset), rownames(params$annotations)),, drop = FALSE]
  }
  if (any(!is.na(params$annotations.genes))) {
    if (any(rownames(data.subset) %notin% rownames(params$annotations.genes))) {
      warning('rownames of input data do not match rownames of annotations.genes, cannot link annotations to data')
    }
    temp.annotations.genes <- params$annotations.genes[match(rownames(data.subset), rownames(params$annotations.genes)),, drop = FALSE]
  }

  if ( all(groupings == FALSE) & all(groupings.genes == FALSE) ) {
    stop('Must specify either "groupings" or "groupings.genes"
       If no groupings are desired consider using myHeatmap instead')
  }

  if (class(groupings) =="character"){
    if ( any(groupings %notin% colnames(temp.annotations))) {
      stop('supplied character vector for groupings not found in sample annotations')}
    factorgroupings <- makefactorgroup(temp.annotations, groupings, specify.gaps = groupings.gaps, return.gaps = TRUE)
    groupings <- factorgroupings$factor.group
    gaps.groupings <- c(factorgroupings$gaps)
    if (any(is.na(groupings))) {
      levs <- levels(groupings)
      groupings <- as.character(groupings); groupings[is.na(groupings)] <- "No_Annot"
      groupings <- factor(groupings, levels = c(levs, "No_Annot"))}
  }else{
    if (class(groupings) == "data.frame") {
      ##order groupings by order of subset
      groupings <- droplevels(groupings[match(colnames(data.subset), rownames(groupings)),1]) %>% as.factor()
      if (any(is.na(groupings))) {
        levs <- levels(groupings)
        groupings <- as.character(groupings); groupings[is.na(groupings)] <- "No_Annot"
        groupings <- factor(groupings, levels = c(levs, "No_Annot"))}
    }
    if (class(groupings) == "factor"){
      if (is.null(names(groupings))==FALSE) {
        groupings <- droplevels(groupings[match(colnames(data.subset), names(groupings))]) %>% as.factor()
        if (any(is.na(groupings))) {
          levs <- levels(groupings)
          groupings <- as.character(groupings); groupings[is.na(groupings)] <- "No_Annot"
          groupings <- factor(groupings, levels = c(levs, "No_Annot"))}
      }else{
        if (length(groupings) != ncol(data.subset)) {
          stop('unnamed factor supplied to groupings not the same length as number of columns in data')
        }
      }  ##do nothing and leave as is
    }
  }



  if (class(groupings.genes) =="character") {
    if (any(groupings.genes %notin% colnames(temp.annotations.genes))) {
      stop('supplied character vector for groupings.genes not found in gene annotations')}
    factorgroupings.genes <- makefactorgroup(temp.annotations.genes, groupings.genes, specify.gaps = groupings.genes.gaps, return.gaps = TRUE)
    groupings.genes <- factorgroupings.genes$factor.group
    gaps.groupings.genes <- c(factorgroupings.genes$gaps)
    if (any(is.na(groupings.genes))) {
      levs <- levels(groupings.genes)
      groupings.genes <- as.character(groupings.genes); groupings.genes[is.na(groupings.genes)] <- "No_Annot"
      groupings.genes <- factor(groupings.genes, levels = c(levs, "No_Annot"))}
  }else{
    if (class(groupings.genes) == "data.frame") {
      ##order groupings by order of subset
      groupings.genes <- droplevels(groupings.genes[match(rownames(data.subset), rownames(groupings.genes)),1]) %>% as.factor()
      if (any(is.na(groupings.genes))) {
        levs <- levels(groupings.genes)
        groupings.genes <- as.character(groupings.genes); groupings.genes[is.na(groupings.genes)] <- "No_Annot"
        groupings.genes <- factor(groupings.genes, levels = c(levs, "No_Annot"))}
    }
    if (class(groupings.genes) == "factor"){
      if (is.null(names(groupings.genes))==FALSE) {
        groupings.genes <- droplevels(groupings.genes[match(rownames(data.subset), names(groupings.genes))]) %>% as.factor()
        if (any(is.na(groupings.genes))) {
          levs <- levels(groupings.genes)
          groupings.genes <- as.character(groupings.genes); groupings.genes[is.na(groupings.genes)] <- "No_Annot"
          groupings.genes <- factor(groupings.genes, levels = c(levs, "No_Annot"))}
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


  if (any(groupings != FALSE, na.rm=T)) {

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

    if (any(groupings.genes == FALSE, na.rm = T)) {  ##if not going on to group genes, see if should be ordered by sample, set clustering of genes
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



  if (any(groupings.genes != FALSE, na.rm=T)) {

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

    if (any(groupings == FALSE, na.rm=T)) {  ##if not grouped by samples, see if should be ordered by gene, set clustering of samples
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
      my_cols <- c(my_cols,params$scale.colors[1],na.color)
    }
    if(is.raw.Ct==FALSE){
      data.subset[which(data.subset1==min(data.subset1))] <- params$scale.range[1]-0.04
      breaks <- c(params$scale.range[1]-0.05,params$scale.range[1]-0.01,breaks)
      my_cols <- c(na.color,params$scale.colors[1],my_cols)
    }
  }


  if (any( rownames(params$annot_samps) %in% colnames(data.subset))) {
    temp.annot_samps <- params$annot_samps
  }else{ temp.annot_samps <- NA}


  if (any( rownames(params$annot_genes) %in% rownames(data.subset))) {
    temp.annot_genes <- params$annot_genes
  }else{ temp.annot_genes <- NA}

  temp.annot_cols <- params$annot_cols

  if (drop.annot.levels == TRUE) {
    if (any(!is.na(temp.annot_samps))) {
      temp.annot_samps[] <- lapply(temp.annot_samps, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_samps <- temp.annot_samps %>% tibble::rownames_to_column("rowN")
      temp.annot_samps <- droplevels(temp.annot_samps[which(temp.annot_samps$rowN %in% colnames(data.subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="rowN")

      spec.cols <- colnames(temp.annot_samps)[colnames(temp.annot_samps) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0 ) {
        for (annot.i in 1:length(spec.cols)) {
          annot <- spec.cols[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_samps[,which(colnames(temp.annot_samps)==annot)])  )]
          if ( sum( levels(temp.annot_samps[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_samps[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_samps[,annot])[levels(temp.annot_samps[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    }


    if (any(!is.na(temp.annot_genes))) {
      temp.annot_genes[] <- lapply(temp.annot_genes, as.factor)
      #subset annot_samps and genes for subset so that annotations will be dropped in heatmap
      temp.annot_genes <- temp.annot_genes %>% tibble::rownames_to_column("rowN")
      temp.annot_genes <- droplevels(temp.annot_genes[which(temp.annot_genes$rowN %in% rownames(data.subset)),]) %>% as.data.frame() %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var="rowN")

      spec.cols <- colnames(temp.annot_genes)[colnames(temp.annot_genes) %in% names(temp.annot_cols)]

      if (length(spec.cols) != 0) {
        for (annot.i in 1:length(spec.cols)) {
          # annot <- colnames(temp.annot_genes)[annot.i]
          annot <- spec.cols[annot.i]
          temp.annot_cols[[which(names(temp.annot_cols)==annot)]] <- temp.annot_cols[[which(names(temp.annot_cols)==annot)]][which(   names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])   %in%   levels(temp.annot_genes[,which(colnames(temp.annot_genes)==annot)])  )]
          if ( sum( levels(temp.annot_genes[,annot]) %in% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])  ) != length(levels(temp.annot_genes[,annot]))) {
            temp.annot_cols[[which(names(temp.annot_cols)==annot)]][c(levels(temp.annot_genes[,annot])[levels(temp.annot_genes[,annot]) %notin% names(temp.annot_cols[[which(names(temp.annot_cols)==annot)]])])] <- "white"
          }
        }
      }
    }

  }



  if (clust.cols == T) {
    tryclustcols <- try(hclust(clust.samps, linkage), silent = T)
    if (class(tryclustcols) == "try-error") {stop('cannot cluster columns, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.cols = F or specify order.by.gene')}
  }

  if (clust.rows == T) {
    tryclustrows <- try(hclust(clust.genes, linkage), silent = T)
    if (class(tryclustrows) == "try-error") {stop('cannot cluster rows, if too many NAs present, set na.fix = T to treat NA values as low expression instead of missing, otherwise set clust.rows = F or specify order.by.sample')}
  }

  if (is.null(main)==TRUE){
    main <- paste("Genes of Interest:",paste(list, collapse = ","))
    main <- paste(main,"\n Method:",method," Linkage:",linkage)}

  pheatmap(data.subset,col=my_cols, breaks=breaks, border_color = border.color, na_col = na.color, clustering_method=linkage,annotation_col=temp.annot_samps, annotation_colors = temp.annot_cols,
           clustering_distance_rows = clust.genes, clustering_distance_cols = clust.samps, main=main,
           cluster_rows = clust.rows, cluster_cols = clust.cols, cutree_rows = row.groups, cutree_cols = col.groups, gaps_row = gaps.row, gaps_col = gaps.col,fontsize_row = fontsize.row, fontsize_col = fontsize.col,
           cellwidth = cell.width, cellheight = cell.height, show_rownames = show.rownames,show_colnames = show.colnames,
           treeheight_row = heightrow ,treeheight_col = heightcol, silent=hide.plot, legend=show.legend, annotation_legend = show.annotations,
           annotation_row=temp.annot_genes, drop_levels = drop.annot.levels, annotation_names_row = annotation.names.row, annotation_names_col = annotation.names.col)

}


#' @export
ExtractMatrix <- function(   ##from clustered heatmap, will extract the exact matrix
  data, ##data heatmap is based off of, (needs to be subset for genes already?)
  heatmap, ##input should be the heatmap saved as a variable, output will be the matrix in the order that the heatmap is output as
  clustered.cols = TRUE, ##assumes that heatmap$tree_col exists
  clustered.rows = TRUE ##assumes that heatmap$tree_row exists
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  attributes <- unlist(lapply(heatmap$gtable$grobs, function(x)(x$name)))
  text.attributes <- grep("GRID.text", attributes)

  ordered <- data

  ##extract columns
  colord <- heatmap$gtable[["grobs"]][[text.attributes[2]]][["label"]]
  ordered <- ordered[,colnames(ordered) %in% colord]
  if(clustered.cols==TRUE){
    ordered <- ordered[,match(colord,colnames(ordered))]
  }

  ##extract rows
  roword <- heatmap$gtable[["grobs"]][[text.attributes[3]]][["label"]]
  ordered <- ordered[rownames(ordered) %in% roword,]
  if(clustered.rows==TRUE){
    ordered <- ordered[match(roword,rownames(ordered)),]
  }

  return(ordered)
}


#' @export
extractClusters <- function(  ##will extract cluster groups from clustered heatmaps
  data, ##should be the same data as what was used to generate the heatmap, used to extract the correct gene order
  heatmap, #this should be the output from pheatmap that was saved
  to.extract = "genes", ##can accept "genes", "samples", or "both"
  nclusters, ## numeric, how many clusters should be extracted. If to.extract is set to "both", numeric vector should be supplied listing number of gene (row clusters) first
  GeneGroup_Name = NULL,
  SampleGroup_Name = NULL,
  all.genes = NULL, #gene names, if annotations should be in placed in larger gene list than matrix in question
  all.samples = NULL #sample names, if annotations should be in placed in larger sample list than matrix in question
){

  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

  if(is.null(all.genes)==TRUE){ all.genes <- rownames(data) }
  if(is.null(all.samples)==TRUE){ all.samples <- colnames(data) }

  if (to.extract == "genes" | to.extract == "both") {
    if (to.extract == "genes") {num_Gene.groups <- nclusters}
    if (to.extract == "both") {num_Gene.groups <- nclusters[1]}
    report.gene.groups <- rep("No_Annot",length(all.genes))
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
  }

  if(to.extract == "samples" | to.extract == "both"){
    if (to.extract == "samples") {num_Sample.groups <- nclusters}
    if (to.extract == "both") {num_Sample.groups <- nclusters[2]}
    report.sample.groups <- rep("No_Annot",length(all.samples))
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
  }

  if(to.extract == "both"){return(c(list(Gene_Groups=report.gene.groups, Sample_Groups=report.sample.groups)))}
  if(to.extract == "genes"){return(report.gene.groups)}
  if(to.extract == "samples"){return(report.sample.groups)}

}




#' @export
extractGaps <- function(  ##similar to extractClusters, will return the positions of gaps created by separating cluster groups
  data, ##should be the same data as what was used to generate the heatmap, used to extract the correct gene order
  heatmap, #this should be the output from pheatmap that was saved
  extractRows = FALSE, ##will extract from the gene clustering
  extractCols = TRUE, ##will extract from the sample clustering
  num_Rows = 5, ##default of 5, should be decided
  num_Cols = 5 ##default of 5, should be decided
){
  if (("matrix" %in% class(data)) != TRUE ) {
    data <- as.matrix(data)
    warning('input data converted to matrix')
  }

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

