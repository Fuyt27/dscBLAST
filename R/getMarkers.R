#' getMarkers
#'
#' @param species
#' options: c('Hs','Mm')
#' @param reference
#' options: c('ecto','meso','endo','exe','pre')
#' @param ref_celltype
#' reference celltype
#' @param ref_dataset
#' reference dataset
#' @param marker_top_n
#' top N markers [default: 20]
#'
#' @return marker.final a list of markers
#' @export
#'
getMarkers <- function(species,reference,ref_celltype,ref_dataset,marker_top_n=20){
  markers_list=loadMarkers(species,reference)
  dataset=names(markers_list)
  markerlist=list()
  for (i in dataset) {
    celltype=unique(markers_list[[i]]$cluster)
    dataset_marker=markers_list[[i]]
    for (j in celltype) {
      celltype_marker=dataset_marker[dataset_marker$cluster==j,]
      markerlist[[i]][[j]]=celltype_marker
    }
  }
  marker.use= markerlist[[ref_dataset]][[ref_celltype]]$gene[1:marker_top_n]
  marker.final =list(marker.use)
  names(marker.final)=c(paste0(ref_dataset,'_',ref_celltype,'_top',marker_top_n))
  return(marker.final)
}

#' getTotalMarkers
#'
#' @param species 
#' options: c('Hs','Mm')
#' @param reference 
#' options: c('ecto','meso','endo','exe','pre')
#' @return markers_list
#' @export

getTotalMarkers <- function(species,reference='all'){
  if(length(reference)==1){if(reference=='all'){reference=c('ecto','meso','endo','exe','pre')}}
  markers_list=list()
  for (i in reference) {
    markers_list[[i]]=loadMarkers(species,i)
  }
  return(markers_list)
  }

#' loadMarkers
#'
#' @param species 
#' options: c('Hs','Mm')
#' @param reference 
#' options: c('ecto','meso','endo','exe','pre')
#' @return markers_list
#' @export
#'
loadMarkers <- function(species,reference){
  if(species=='Hs'){
    if (reference=='ecto') {
      markers_list=markers_list_ectoH
    }
    else if (reference=='endo') {
      markers_list=markers_list_endoH
    }
    else if (reference=='meso') {
      markers_list=markers_list_mesoH
    }
    else if (reference=='pre') {
      markers_list=markers_list_preH
    }
    else if (reference=='exe') {
      markers_list=markers_list_exeH
    }
  }
  else if (species=='Mm'){
    if (reference=='ecto') {
      markers_list=markers_list_ectoM
    }
    else if (reference=='endo') {
      markers_list=markers_list_endoM
    }
    else if (reference=='meso') {
      markers_list=markers_list_mesoM
    }
    else if (reference=='pre') {
      markers_list=markers_list_preM
    }
    else if (reference=='exe') {
      markers_list=markers_list_exeM
    }
  }
  return(markers_list)
}