####Parameter setup and Updating####
#' @export
parameters <- list(scale.range = c(-1,1),
                   scale.colors = c("blue","black","yellow"),
                   n.colors.range = 100,
                   annotations = NA,
                   annot_samps = NA,
                   annotations.genes = NA,
                   annot_genes = NA,
                   annot_cols = NA,
                   expression_gradient.colors = c("blue","lightblue","gray","indianred","firebrick"))




#' @export
initiate_params <- function(params = parameters){
  assign("params",params,envir =  .GlobalEnv)
  message('params list object containing parameters and annotations has been set to default values')
}



#' @export
set_scale.range <- function(range){
  if(length(range) == 2 & class(range) == "numeric" & range[2] > range[1]){
    params$scale.range <<- range}
  else{stop('range must be a numeric of length 2 and the second element must be greater than the first')}}



#' @export
set_scale.colors <- function(colors){
  params$scale.colors <<- colors}



#' @export
set_n.colors.range <- function(n){params$n.colors.range <<- n}



#' @export
set_annotations <- function(annotations){
  if (any(is.na(annotations))) {
    if ( (length(annotations) == 1) & (is.null(nrow(annotations)) == TRUE) ) {
      params$annotations <<- NA
    }else{
      if (length(.row_names_info(annotations, type = 0)) != nrow(annotations) ) {
        stop('No rownames present in provided annotations. Rownames of annotations are linked to colnames of data and must be provided for proper integration')}
      params$annotations <<- annotations
      params$annotations[is.na(params$annotations)] <<- "No_Annot"
      warning('NAs present in annotations, converted to "No_Annot" for functionality purposes')
    }
  }
  params$annotations <<- annotations
}



#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr full_join
#' @export
update_annotations <- function(   ###should make an option to full join if stuff thats being added doesn't cover everything....
  annotation,
  values
){
  temp.annotations <- params$annotations

  if (any(annotation %in% colnames(temp.annotations))) {   ##if being updated, removes prior instance and will be added back in below
    temp.annotations <- temp.annotations[,-which(colnames(temp.annotations) %in% annotation), drop = FALSE]}


  if (is.null(dim(values)) == TRUE) {  ##vector or factor, not matrix or dataframe
    if (is.null(names(values)) == TRUE) {  ##vector with no names supplied, assumes order and length is the same
      if (length(values) != nrow(temp.annotations)) {stop('unnamed vector supplied for values is not the same length as the number of rows in annotations, please set names to sample names for proper integration ')
      }else{
        temp.annotations$V1 <- values
        colnames(temp.annotations)["V1"] <- annotation
      }
      new.cols <- data.frame(V1 = values, RowNames = names(values)); colnames(new.cols)["V1"] <- annotation
      temp.annotations <- temp.annotations %>% tibble::rownames_to_column("RowNames")
      suppressMessages(temp.annotations <- dplyr::full_join(temp.annotations, new.cols) %>% tibble::column_to_rownames("RowNames"))
    }
  }else{
    if (length(annotation) != ncol(values)) {stop('length of supplied annotation names not equal to number of value column provided')}

    if (length(.row_names_info(values, type = 0)) != nrow(values) ) {
      stop('No rownames present in provided values. Rownames must be provided to properly update annotations')}


    new.cols <- as.data.frame(values); colnames(new.cols) <- annotation; new.cols <- new.cols %>% tibble::rownames_to_column("RowNames")
    temp.annotations <- temp.annotations %>% tibble::rownames_to_column("RowNames")
    suppressMessages(temp.annotations <- dplyr::full_join(temp.annotations, new.cols) %>% tibble::column_to_rownames("RowNames"))
  }

  if (any(is.na(temp.annotations))) {
    for (i in 1:ncol(temp.annotations)) {
      if (any(is.na(temp.annotations[,i]))) {
        if (class(temp.annotations[,i]) == "factor") {
          levs <- levels(temp.annotations[,i])
          temp.annotations[,i] <- as.character(temp.annotations[,i])
          temp.annotations[,i][is.na(temp.annotations[,i])] <- "No_Annot"
          warning('NAs present in updated annotations, converted to "No_Annot" for functionality purposes')
          temp.annotations[,i] <- factor(temp.annotations[,i], levels =c(levs,"No_Annot"))
        }else{
          temp.annotations[,i][is.na(temp.annotations[,i])] <- "No_Annot"
          warning('NAs present in updated annotations, converted to "No_Annot" for functionality purposes')
        }
      }
    }
  }

  params$annotations <<- temp.annotations

}




#' @export
set_annot_samps <- function(annotations= NULL){
  if (is.null(annotations) == TRUE) { params$annot_samps <<- params$annotations
  }else{
    if (all(is.na(annotations)))  { params$annot_samps <<- NA
    }else{
      if (sum(annotations %in% colnames(params$annotations)) != length(annotations)) {
        stop('one or more of the supplied list of annotations cannot be found in annotations')
      }
      params$annot_samps <<- params$annotations[,which(colnames(params$annotations) %in% annotations), drop = FALSE]
      params$annot_samps <<- params$annot_samps[,match(annotations, colnames(params$annot_samps)), drop = FALSE]
    }
  }
}



#' @export
set_annotations.genes <- function(annotations.genes){

  if (any(is.na(annotations.genes))) {
    if ( (length(annotations.genes) == 1) & (is.null(nrow(annotations.genes)) == TRUE) ) {
      params$annotations.genes <<- NA
    }else {
      if (length(.row_names_info(annotations.genes, type = 0)) != nrow(annotations.genes) ) {
        stop('No rownames present in provided annotations. Rownames of annotations are linked to rownames of data and must be provided for proper integration')}

      params$annotations.genes <<- annotations.genes
      params$annotations.genes[is.na(params$annotations.genes)] <<- "No_Annot"

      warning('NAs present in annotations.genes, converted to "No_Annot" for functionality purposes')
    }
  }
  params$annotations.genes <<- annotations.genes
}


#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr full_join
#' @export
update_annotations.genes <- function(   ###should make an option to full join if stuff thats being added doesn't cover everything....
  annotation,
  values
){
  temp.annotations.genes <- params$annotations.genes

  if (any(annotation %in% colnames(temp.annotations.genes))) {   ##if being updated, removes prior instance and will be added back in below
    temp.annotations.genes <- temp.annotations.genes[,-which(colnames(temp.annotations.genes) %in% annotation), drop = FALSE]}


  if (is.null(dim(values)) == TRUE) {  ##vector or factor, not matrix or dataframe
    if (is.null(names(values)) == TRUE) {  ##vector with no names supplied, assumes order and length is the same
      if (length(values) != nrow(temp.annotations.genes)) {stop('unnamed vector supplied for values is not the same length as the number of rows in annotations, please set names to sample names for proper integration ')
      }else{
        temp.annotations.genes$V1 <- values
        colnames(temp.annotations.genes)["V1"] <- annotation
      }
      new.cols <- data.frame(V1 = values, RowNames = names(values)); colnames(new.cols)["V1"] <- annotation
      temp.annotations.genes <- temp.annotations.genes %>% tibble::rownames_to_column("RowNames")
      suppressMessages(temp.annotations.genes <- dplyr::full_join(temp.annotations.genes, new.cols) %>% tibble::column_to_rownames("RowNames"))
    }
  }else{
    if (length(annotation) != ncol(values)) {stop('length of supplied annotation names not equal to number of value column provided')}

    if (length(.row_names_info(values, type = 0)) != nrow(values) ) {
      stop('No rownames present in provided values. Rownames must be provided to properly update annotations')}


    new.cols <- as.data.frame(values); colnames(new.cols) <- annotation; new.cols <- new.cols %>% tibble::rownames_to_column("RowNames")
    temp.annotations.genes <- temp.annotations.genes %>% tibble::rownames_to_column("RowNames")
    suppressMessages(temp.annotations.genes <- dplyr::full_join(temp.annotations.genes, new.cols) %>% tibble::column_to_rownames("RowNames"))
  }


  if (any(is.na(temp.annotations.genes))) {
    for (i in 1:ncol(temp.annotations.genes)) {
      if (any(is.na(temp.annotations.genes[,i]))) {
        if (class(temp.annotations.genes[,i]) == "factor") {
          levs <- levels(temp.annotations.genes[,i])
          temp.annotations.genes[,i] <- as.character(temp.annotations.genes[,i])
          temp.annotations.genes[,i][is.na(temp.annotations.genes[,i])] <- "No_Annot"
          warning('NAs present in updated annotations, converted to "No_Annot" for functionality purposes')
          temp.annotations.genes[,i] <- factor(temp.annotations.genes[,i], levels =c(levs,"No_Annot"))
        }else{
          temp.annotations.genes[,i][is.na(temp.annotations.genes[,i])] <- "No_Annot"
          warning('NAs present in updated annotations, converted to "No_Annot" for functionality purposes')
        }
      }
    }
  }

  params$annotations.genes <<- temp.annotations.genes

}




#' @export
set_annot_genes <- function(annotations = NULL){
  if (is.null(annotations) == TRUE) { params$annot_genes <<- params$annotations.genes
  }else{
    if (all(is.na(annotations))) { (params$annot_genes <<- NA)
    }else{
      if (sum(annotations %in% colnames(params$annotations.genes)) != length(annotations)) {
        stop('one or more of the supplied list of annotations cannot be found in annotations.genes')
      }
      params$annot_genes <<- params$annotations.genes[,which(colnames(params$annotations.genes) %in% annotations),drop=FALSE]
      params$annot_genes <<- params$annot_genes[,match(annotations, colnames(params$annot_genes)), drop = FALSE]
    }

  }

}




#' @export
set_annot_cols <- function(annot_cols){
  params$annot_cols <<- annot_cols
  if (all(is.na(annot_cols))) { params$annot_cols <<- NA}
}




#' @export
update_annot_cols <- function(  ##right now must be one at a time
  annotation,
  values.list
){
  temp.annot_cols <- params$annot_cols

  if (annotation %in% names(temp.annot_cols)) {   ##if being updated, removes prior instance and will be added back in below
    temp.annot_cols <- temp.annot_cols[-which(names(temp.annot_cols) %in% annotation)]}

  if (all(!is.na(values.list))) {
    temp.annot_cols$V1 <- values.list
    names(temp.annot_cols)[names(temp.annot_cols)=="V1"] <- annotation
  }

  params$annot_cols <<- temp.annot_cols
}


#' @importFrom RColorBrewer brewer.pal
#' @export
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
