#' Subset data
#'
#' @param auc
#' correlation bewtween query and ref celltypes(results from RUN_dscBLAST)
#' @param top_n
#' top n correlated cell_type
#' default: 3
#'
#' @param ref_dir
#' directory containing count and meta
#' @return ref_list a list of reference single cell experiment objects
#' @export
#'
subset_data <- function(auc,top_n,ref_dir){
  aurocs2 = auc %>% as.data.frame()
  label = colnames(aurocs2)
  item = list()
  if(dim(aurocs2)[2]==1){aurocs2$V2=0}
  for (i in label){
    aurocs2 = aurocs2[order(-aurocs2[[i]]),]
    item[[i]] = rownames(aurocs2)[1:top_n]
  }
  select = unlist(item)

  dataset =unique(colsplit(select,'_',names = c('c1','c2'))$c1)
  select_sample_cell =colsplit(select,'_',names = c('c1','c2'))$c2
  ref_list=list()
  for (i in dataset) {
    if(i%in%c('ectoM','mesoM','endoM','preM','exeM','ectoH','mesoH','endoH','preH','exeH')){
      dge <- readRDS(paste0(ref_dir,'/','data_',i,'.rds'))
      meta <-readRDS(paste0(ref_dir,'/','meta_',i,'.rds'))
      ref_list[[i]] <- CreateSeuratObject(dge);
      ref_list[[i]]@meta.data <- meta

      ref_list[[i]] <- subset(ref_list[[i]],sample_cell %in% select_sample_cell)
      ref_list[[i]] = as.SingleCellExperiment(ref_list[[i]])
    }
  }
  return(ref_list)
}

#' Subset marker
#'
#' @param ref_list
#' reference list
#' @param marker_top_n
#' top n markers
#' default: 20
#'
#' @return ref_marker_list
#' @export
#'
subset_marker <- function(ref_list,marker_top_n=20){
  ref_marker_list=list()
  for (i in names(ref_list)) {
    ref_marker_list[[i]] <- get_topN_markers(ref_list,i,marker_top_n)  }
  return(ref_marker_list)
}

#' Get topN markers from each dataset
#'
#' @param ref_list
#' ref data list
#' @param reference
#' 'ectoM','mesoM','endoM','preM','exeM','ectoH','mesoH','endoH','preH','exeH'
#' @param marker_top_n
#' default: 20
#'
#' @return topN markers from chosen dataset
#' @export
#' @importFrom dplyr %>% slice_max group_by
get_topN_markers <- function(ref_list,reference,marker_top_n){
  if(reference=='ectoM'){
    markers_list = markers_list_ectoM[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]}
  if(reference=='endoM'){
    markers_list = markers_list_endoM[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='mesoM'){
    markers_list = markers_list_mesoM[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='exeM'){
    markers_list = markers_list_exeM[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='preM'){
    markers_list = markers_list_preM[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='ectoH'){
    markers_list = markers_list_ectoH[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]}
  if(reference=='endoH'){
    markers_list = markers_list_endoH[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='mesoH'){
    markers_list = markers_list_mesoH[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='exeH'){
    markers_list = markers_list_exeH[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  if(reference=='preH'){
    markers_list = markers_list_preH[unique(ref_list[[reference]]$dataset)]
    total_list = list()
    for (j in 1:length(markers_list)) {
      set=markers_list[[j]]  %>% as.data.frame()
      set=set[set$cluster%in%unique(ref_list[[reference]]$cell_type),]
      topN <- set%>% group_by(cluster) %>% slice_max(avg_log2FC,n=marker_top_n)
      total_list[[j]]=topN$gene
    }
    marker_topN = unlist(total_list)
    marker_topN = marker_topN[!duplicated(marker_topN)]
    marker_topN = marker_topN[!is.na(marker_topN)]
  }
  return(marker_topN)
}

#' subRight character
#'
#' @param x
#' a character
#' @param n
#' order from right[default: 1]
#' @return sub_right one character
#' @export
#'
substrRight <- function(x, n=1){
  sub_right=substr(x, nchar(x)-n+1, nchar(x))
  return(sub_right)
}

#' Calculate auroc with stage info
#'
#' @param sce
#' a single cell experiment object
#' @param pretrained_model
#' a pretrained_model list
#' @param query_species
#' options: c('Hs','Mm')
#'
#' @return correlation bewtween query and ref celltypes
#' @export
#'
auc_calculate_withstage <- function(sce,pretrained_model,query_species,highlight){
  auc_list=list()
  if (!is.list(sce)) {
    for (i in names(pretrained_model)) {
      if (dim(pretrained_model[[i]])[2]==1){next}
      auc_list[[i]] = MetaNeighborUS(
        trained_model = pretrained_model[[i]], dat = sce,
        study_id = sce$study_id, cell_type = sce$cell_type,
        fast_version = TRUE,node_degree_normalization=F)
      colnames(auc_list[[i]])=paste0(i,'_',colnames(auc_list[[i]]))
      auc_list[[i]]=t(auc_list[[i]])
    }
    target.num=length(unique(sce$cell_type))
    total=data.frame(matrix(nrow=0,ncol =target.num))
    for (i in 1:length(auc_list)) {
      tmp=auc_list[[i]]
      total=rbind(total,tmp)
    }
    if(highlight){
      highlight_celltype=colnames(auc[['auc_highlight']])
        if(length(highlight_celltype)>1){
          highlight=total[,highlight_celltype]}
        else{auc_tmp=matrix(nrow = dim(total)[1],ncol = 1);colnames(auc_tmp)=highlight_celltype;rownames(auc_tmp)=rownames(total);auc_tmp[,highlight_celltype]=total[,highlight_celltype];highlight=auc_tmp}
        total_auc =list(auc_list=auc_list,auc_total=total,auc_highlight=highlight)
        auc_final=list(auc_list=auc_list,auc_total=total,auc_highlight=highlight)
        }
    else{auc_final=list(auc_list=auc_list,auc_total=total)}
  }
  else{
    for (i in names(pretrained_model)) {
      if (dim(pretrained_model[[i]])[2]==1){next}
      use_species=ifelse(substrRight(i,1) == "H", 'Hs', 'Mm')
      if (use_species==query_species) {
        auc_list[[i]] = MetaNeighborUS(
          trained_model = pretrained_model[[i]], dat = sce[[query_species]],
          study_id = sce[[query_species]]$study_id, cell_type = sce[[query_species]]$cell_type,
          fast_version = TRUE,node_degree_normalization=F)
        colnames(auc_list[[i]])=paste0(i,'_',colnames(auc_list[[i]]))
        auc_list[[i]]=t(auc_list[[i]])}
      else{
        trans_species=ifelse(query_species == "Hs", 'Mm', 'Hs')
        auc_list[[i]] = MetaNeighborUS(
          trained_model = pretrained_model[[i]], dat = sce[[trans_species]],
          study_id = sce[[trans_species]]$study_id, cell_type = sce[[trans_species]]$cell_type,
          fast_version = TRUE,node_degree_normalization=F)
        colnames(auc_list[[i]])=paste0(i,'_',colnames(auc_list[[i]]))
        auc_list[[i]]=t(auc_list[[i]])
      }
    }
    target.num=length(unique(sce[[query_species]]$cell_type))
    total=data.frame(matrix(nrow=0,ncol =target.num))
    for (i in 1:length(auc_list)) {
      tmp=auc_list[[i]]
      total=rbind(total,tmp)
    }
    if(highlight){
      highlight_celltype=colnames(auc[['auc_highlight']])
      if(length(highlight_celltype)>1){
        highlight=total[,highlight_celltype]}
      else{auc_tmp=matrix(nrow = dim(total)[1],ncol = 1);colnames(auc_tmp)=highlight_celltype;rownames(auc_tmp)=rownames(total);auc_tmp[,highlight_celltype]=total[,highlight_celltype];highlight=auc_tmp}
      total_auc =list(auc_list=auc_list,auc_total=total,auc_highlight=highlight)
      auc_final=list(auc_list=auc_list,auc_total=total,auc_highlight=highlight)
    }
    else{auc_final=list(auc_list=auc_list,auc_total=total)}


  }
  return(auc_final)
}

#' Gene transfer
#'
#' @param expression_matrix
#' a gene expression matrix
#' @param query_species
#' Hs/Mm
#' @return a transfered gene expression matrix(human to mouse or mouse to human)
#' @export
#'

transfer_gene <- function(expression_matrix,query_species){
  if(query_species=='Hs'){
    Mouse_gene <- gene_transfer[gene_transfer$Human_gene%in%rownames(expression_matrix),]
    Mouse_gene<-Mouse_gene[!duplicated(Mouse_gene$Mouse_gene),]
    Mouse_gene<-Mouse_gene[!Mouse_gene$Mouse_gene=='',]
    expression_matrix <- expression_matrix[Mouse_gene$Human_gene,]
    rownames(expression_matrix) <-Mouse_gene$Mouse_gene
    return(expression_matrix)
  }
  if(query_species=='Mm'){
    Human_gene <- gene_transfer[gene_transfer$Mouse_gene%in%rownames(expression_matrix),]
    Human_gene<-Human_gene[!duplicated(Human_gene$Human_gene),]
    Human_gene<-Human_gene[!Human_gene$Human_gene=='',]
    expression_matrix <- expression_matrix[Human_gene$Mouse_gene,]
    rownames(expression_matrix) <-Human_gene$Human_gene
    return(expression_matrix)
  }
}

#' Cell qc
#'
#' @param expression_matrix
#' a gene-cell expression matrix
#' @param query_species
#' options: c('Hs','Mm')
#' @param metadata
#' meta data 
#' @param cell_type
#' the cell type info of query cells
#' @param downsample
#' set number of cells
#' default: 50000
#' @return a seurat object of cells passing qulity control

#' @export
#' @importFrom Seurat CreateSeuratObject PercentageFeatureSet

Cell_qc<-function(expression_matrix,query_species,metadata,cell_type,downsample) {
  if (dim(expression_matrix)[1]==0 ||dim(expression_matrix)[2]==0) {
    stop("`expression_matrix`should be a gene-cell matrix")
  }
  seob=CreateSeuratObject(counts =expression_matrix, min.cells = 3, min.features = 200)
  seob$orig.ident='SeuratObject'
  seob@meta.data=metadata
  seob$cell_type=cell_type
  message('--step1:detect MT genes--')
  if (query_species=='Mm') {
    seob$percent.mt=PercentageFeatureSet(seob,pattern = '^mt')
  }
  else if(query_species=="Hs"){
    seob$percent.mt=PercentageFeatureSet(seob,pattern = '^MT')
  }
  seob<-subset(seob,percent.mt <25)
  message('--step2:downsample--')
  if(dim(seob)[2]<=downsample){
    return(seob)
  }
  else{
    set.seed(1)
    seob<-seob[,sample(dim(seob)[2],downsample)]
    return(seob)}
}

#' Create matrix for big dataset
#'
#' @param mat
#' a sparse matrix
#' @return a matrix
#' @export
#'
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x

  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}


#' train model with stage info
#'
#' @param var_genes
#' variable genes
#' @param dat
#' a single cell experiment object
#'
#' @param i
#' a number[default:1]
#' @param study_id
#' ref dataset info
#'
#' @param cell_type
#' ref celltype info
#'
#' @param stage_id
#' ref stage info
#'
#' @return a pretrained model with stage info
#'
#' @export
#'

trainModel_stage<-function (var_genes, dat, i = 1, study_id, cell_type,stage_id)
{
  dat <- SummarizedExperiment::assay(dat, i = i)
  samples <- colnames(dat)
  if (length(study_id) != length(samples)) {
    stop("study_id length does not match number of samples")
  }
  if (length(cell_type) != length(samples)) {
    stop("cell_type length does not match number of samples")
  }
  matching_vargenes <- match(rownames(dat), var_genes)
  matching_vargenes_count <- sum(!is.na(matching_vargenes))
  if (matching_vargenes_count < 2) {
    stop("matching_vargenes should have more than 1 matching genes!",
         call. = TRUE)
  }
  else if (matching_vargenes_count < 5) {
    warning("matching_vargenes should have more matching genes!",
            immediate. = TRUE)
  }
  dat <- dat[!is.na(matching_vargenes), ]
  study_id <- as.character(study_id)
  cell_type <- as.character(cell_type)
  dat <- normalize_cols(dat)
  label_matrix <- design_matrix(paste(study_id, cell_type, stage_id,
                                      sep = "|"))
  result <- dat %*% label_matrix
  result <- rbind(n_cells = colSums(label_matrix), result)
  return(result)
}



#' Scale a matrix
#'
#' @param M
#' a matrix
#'
#' @return a scaled matrix
#' @export
#'
scale_cols <- function(M) {
  cm <- colMeans(M)
  cnorm <- 1 / sqrt(colSums(M**2) - nrow(M) * cm**2)
  matrixStats::t_tx_OP_y(matrixStats::t_tx_OP_y(M, cm, "-"), cnorm, "*")
}


#' Normalize matrix
#'
#' @param M a matrix
#' @param ranked
#' default: T
#'
#' @return a normalized matrix
#' @export
#'
normalize_cols <- function(M, ranked = TRUE) {
  result <- as_matrix(M)
  if (ranked) {
    result <- matrixStats::colRanks(result, ties.method = "average",
                                    preserveShape = TRUE)
  }
  result <- scale_cols(result)
  dimnames(result) <- dimnames(M)
  return(result)
}


#' Transform a vector with cell_type labels into a binary matrix
#'
#' @param cell_type
#' ref celltype info
#'
#' @return a binary matix added celltype info
#' @export
#'
design_matrix <- function(cell_type) {
  cell_type <- as.factor(cell_type)
  if (length(levels(cell_type)) > 1) {
    result <- stats::model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- levels(cell_type)
  return(result)
}


