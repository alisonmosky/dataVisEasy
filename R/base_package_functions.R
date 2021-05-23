####Standard Functions, Functions Used within the Package, Data Assessment and Manipulation####

#' @export
'%notin%' <- function(x, table){ !(match(x, table, nomatch = 0) > 0)}




#' @export
assessScale <- function(     ###Returns percent of data above, below, and within range of scale for heatmap
  data ##data to be assessed,
){
  above <- sum(data > params$scale.range[2], na.rm=T)/sum(!is.na(data))
  within <- sum(data <= params$scale.range[2] & data >= params$scale.range[1], na.rm=T)/sum(!is.na(data))
  below <- sum(data < params$scale.range[1], na.rm=T)/sum(!is.na(data))
  report <- c(below, within, above); names(report) <- c("Percent Below Range", "Percent Within Range","Percent Above Range")
  return(report)
}




#' @export
subsetGenes <- function(   ##simply pulls list of genes out of dataset into a new matrix
  data,  ##this will be a matrix of values, with genes in the rows and samples in the columns
  list,   ##list of gene or genes that will be extracted from the dataset, exact matches
  order.by = NULL, ##gene to order by
  exact = TRUE ##if true, it'll find exact matches, if false, will use grep and find anything with the term in it
){

  if(exact==TRUE){
    subset <- data[which(rownames(data) %in% list),];
    if (nrow(subset) == 0 ) {stop('exact matches for list not found in rownames data')}
  }

  if (exact==FALSE){
    subset <- data[grep(paste(list, collapse = "|"),rownames(data)),]
    if (nrow(subset) == 0 ) {stop('inexact matches for list not found in rownames data')}
  }

  if (is.null(order.by)==FALSE){
    subset <- subset[,order(data[which(rownames(data) %in% order.by),],na.last = F)]
  }
  return(subset)
}




#' @export
subsetSamples <- function(  ##subset samples out of matrix based on metadata group
  data, ##matrix of values, with genes in the rows and samples in the columns
  group, ##group from which to subset, vector that is same length as columns
  take.out ##which part of the group to extract, one or more of the factors in the "group"
){

  temp.annotations <- params$annotations

  if (group %in% colnames(temp.annotations)) {

    # if (sum(colnames(data) %notin% rownames(temp.annotations)) != 0 ) {
    if (any(colnames(data) %notin% rownames(temp.annotations))) {
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




#' @importFrom grDevices col2rgb colorRamp colorRampPalette rgb
#' @export
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


#' @export
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
