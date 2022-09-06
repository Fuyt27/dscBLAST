#' Title: RUN_dscBLAST.default
#' @param sce
#' a single cell experiment object
#' @param query_species
#' options: c('Hs','Mm')
#' @param reference
#' options: c('ecto','meso','endo','exe','pre')
#' @return aurocs
#' correlation bewtween query and ref celltypes
#' @export
#'
#' @importFrom MetaNeighbor MetaNeighborUS

RUN_dscBLAST.default<- function(sce,query_species,reference) {
  if (query_species=='Mm') {
    if (reference=='ecto') {
      aurocs = MetaNeighborUS(
        trained_model =m.ecto.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('ectoM_',colnames(aurocs))
      }
    else if(reference=='meso'){
      aurocs = MetaNeighborUS(
        trained_model =m.meso.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('mesoM_',colnames(aurocs))
      }
    else if(reference=='endo'){
      aurocs = MetaNeighborUS(
        trained_model =m.endo.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('endoM_',colnames(aurocs))
    }
    else if(reference=='exe'){
      aurocs = MetaNeighborUS(
        trained_model =m.exe.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('exeM_',colnames(aurocs))
    }
    else if(reference=='pre'){
      aurocs = MetaNeighborUS(
        trained_model =m.pre.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('preM_',colnames(aurocs))
    }
    }
  else if(query_species=="Hs"){
    if (reference=='ecto') {
      aurocs = MetaNeighborUS(
        trained_model =h.ecto.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('ectoH_',colnames(aurocs))
      }
    else if(reference=='meso'){
      aurocs = MetaNeighborUS(
        trained_model =h.meso.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('mesoH_',colnames(aurocs))
      }
    else if(reference=='endo'){
      aurocs = MetaNeighborUS(
        trained_model =h.endo.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('endoH_',colnames(aurocs))
    }
    else if(reference=='exe'){
      aurocs = MetaNeighborUS(
        trained_model =h.exe.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('exeH_',colnames(aurocs))
    }
    else if(reference=='pre'){
      aurocs = MetaNeighborUS(
        trained_model =h.pre.pretrained_model, dat = sce,
        study_id =sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE)
      colnames(aurocs)=paste0('preH_',colnames(aurocs))
    }
  }
  return(t(aurocs))
}

#' Title:RUN dscBLAST
#'
#' @param sce
#' a single cell experiment object
#' @param query_species
#' options: c('Hs','Mm')

#' @param reference
#' options: c('all','ecto','meso','endo','exe','pre')
#' default: 'all'
#'
#' @param highlight_celltype
#' celltype to highlight
#'
#' @return auroc
#' correlation bewtween query and ref celltypes
#' @export
#'
RUN_dscBLAST<- function(sce,query_species,reference='all',highlight_celltype=NULL) {
if (!is.list(sce)) {
    message('--Calculating AUC--')
    if(length(reference)==1){
    if(reference=='all'){reference=c('ecto','meso','endo','exe','pre')}}
    aurocs =list()
    for (i in reference){
    aurocs[[i]]=RUN_dscBLAST.default(sce,query_species,i)
    }
    names(aurocs)=reference
    target.num=length(unique(sce$cell_type))
    total=data.frame(matrix(nrow=0,ncol =target.num))
    for (i in reference) {
      tmp=aurocs[[i]]
      total=rbind(total,tmp)
    }
    if(!is.null(highlight_celltype)){
      used_celltype=paste0('Query|',highlight_celltype)
      if(length(highlight_celltype)>1){
        highlight=total[,used_celltype]}
      else{auc_tmp=matrix(nrow = dim(total)[1],ncol = 1);colnames(auc_tmp)=used_celltype;rownames(auc_tmp)=rownames(total);auc_tmp[,used_celltype]=total[,used_celltype];highlight=auc_tmp}
      total_auc =list(auc_list=aurocs,auc_total=total,auc_highlight=highlight)
      }
    else{total_auc =list(auc_list=aurocs,auc_total=total)}
 }
  else{
    if(length(reference)==1){
      if(reference=='all'){reference=c('ecto','meso','endo','exe','pre')}}
      message('--Calculating AUC--')
      aurocs1 =list()
      for (i in reference){
        aurocs1[[i]]=RUN_dscBLAST.default(sce[[query_species]],query_species,i)}
      aurocs2 =list()
      for (i in reference){
        trans_species=ifelse(query_species == "Hs", 'Mm', 'Hs')
        aurocs2[[i]]=RUN_dscBLAST.default(sce[[trans_species]],trans_species,i)
      }
      names(aurocs1)=reference
      names(aurocs2)=reference
      aurocs<-list(aurocs1,aurocs2)
      names(aurocs)=c(query_species,trans_species)

      target.num=length(unique(sce[[query_species]]$cell_type))

      total_1=data.frame(matrix(nrow=0,ncol =target.num))
      for (i in 1:length(aurocs1)) {
        tmp=aurocs1[[i]]
        total_1=rbind(total_1,tmp)
      }

      total_2=data.frame(matrix(nrow=0,ncol =target.num))
      for (i in 1:length(aurocs2)) {
        tmp=aurocs2[[i]]
        total_2=rbind(total_2,tmp)
      }
      total=rbind(total_1,total_2)
      if(!is.null(highlight_celltype)){
        used_celltype=paste0('Query|',highlight_celltype)
        if(length(highlight_celltype)>1){
          highlight=total[,used_celltype]}
        else{auc_tmp=matrix(nrow = dim(total)[1],ncol = 1);colnames(auc_tmp)=used_celltype;rownames(auc_tmp)=rownames(total);auc_tmp[,used_celltype]=total[,used_celltype];highlight=auc_tmp}
        total_auc =list(auc_list=aurocs,auc_total=total,auc_highlight=highlight)
      }else{total_auc =list(auc_list=aurocs,auc_total=total)}
    }

  message('--finished--')

  return(total_auc)
}


