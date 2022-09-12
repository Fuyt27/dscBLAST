#' Title:create dscBLAST object
#'
#' @param expression_profile
#' cell-gene matrix
#' @param query_species
#' options: c('Hs','Mm')
#' @param ref_species
#' options: c('both','single')
#' default: both
#'
#' @param metadata
#' meta data containing a `cell-type` column
#'
#' @param batch
#' batch info in metadata needed in SCT transform
#' if you offer a unnormalized matrix and want to normalize it by batch
#' in this case,the metadata is expected to contain a `batch` column
#'
#' @param downsample
#' set number of cells
#' default: 50000
#'
#' @param mtx.type
#' set the matrix type
#' if you offer a normalized matrix,choose 'normalized' option
#' default: raw
#'
#' @return sce a single cell experiment object or a list of single cell experiment  objects
#'
#'
#' @export
#'
#' @importFrom Seurat PercentageFeatureSet as.SingleCellExperiment CreateSeuratObject SCTransform SplitObject
create_dscBLASTobject <- function(expression_profile,query_species,ref_species='both',metadata,cell_type,batch='default',downsample=50000,mtx.type='raw') {
  if(ref_species!='both'){
    seob=Cell_qc(expression_profile,query_species,metadata,cell_type,downsample)
    if(mtx.type=='raw'){
    message('--step3:SCT transform--')
    if(batch=='default'){
      seob<- SCTransform(
        seob,
        variable.features.n = 3000,
        vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
        verbose = T)
    }
    else{
    seob_list = SplitObject(seob,split.by = batch)
    message('--step3:SCT transform--')
    for(i in 1:length(seob_list)){
      seob_list[[i]] <- SCTransform(
        seob_list[[i]],
        variable.features.n = 3000,
        vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
        verbose = T)
    }
    seob = merge(seob_list[[1]],seob_list[-1])}
    message('--step4:creating dscBLAST object--')
    seob_sct = CreateSeuratObject(counts = seob@assays[["SCT"]]@data)
    seob_sct@meta.data=seob@meta.data
    sce=as.SingleCellExperiment(seob_sct)
    sce$study_id='Query'
    return(sce)
    message('--finished--')}
    else if (mtx.type=='normalized'){
      message('--step3:creating dscBLAST object--')
      seob_prenorm= CreateSeuratObject(counts = seob@assays[["RNA"]]@data)
      seob_prenorm@meta.data=seob@meta.data
      sce=as.SingleCellExperiment(seob_prenorm)
      sce$study_id='Query'
      return(sce)
    }
    }
  else{
    seob=Cell_qc(expression_profile,query_species,metadata,cell_type,downsample)
    if(mtx.type=='raw'){
    if(batch=='default'){
      message('--step3:SCT transform--')
      seob<- SCTransform(
        seob,
        variable.features.n = 3000,
        vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
        verbose = T)
    }
    else{
      seob_list = SplitObject(seob,split.by = batch)
      message('--step3:SCT transform--')
      for(i in 1:length(seob_list)){
        seob_list[[i]] <- SCTransform(
          seob_list[[i]],
          variable.features.n = 3000,
          vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
          verbose = T)
      }
      seob = merge(seob_list[[1]],seob_list[-1])}
    message('--step4:creating dscBLAST object--')
    counts = seob@assays[["SCT"]]@data
    seob_sct = CreateSeuratObject(counts = counts)
    seob_sct@meta.data=seob@meta.data
    sce1=as.SingleCellExperiment(seob_sct)
    sce1$study_id='Query'

    counts2<-transfer_gene(counts,query_species)
    seob_sct = CreateSeuratObject(counts = counts2)
    seob_sct@meta.data=seob@meta.data
    sce2=as.SingleCellExperiment(seob_sct)
    sce2$study_id='Query'
    sce_list=list(sce1,sce2);names(sce_list)=c(query_species,ifelse(query_species == "Hs", 'Mm', 'Hs'))
    return(sce_list)}
    else if(mtx.type=='normalized'){
      message('--step3:creating dscBLAST object--')
      counts = seob@assays[["RNA"]]@data
      seob_prenorm = CreateSeuratObject(counts = counts)
      seob_prenorm@meta.data=seob@meta.data
      sce1=as.SingleCellExperiment(seob_prenorm)
      sce1$study_id='Query'

      counts2<-transfer_gene(counts,query_species)
      seob_prenorm = CreateSeuratObject(counts = counts2)
      seob_prenorm@meta.data=seob@meta.data
      sce2=as.SingleCellExperiment(seob_prenorm)
      sce2$study_id='Query'
      sce_list=list(sce1,sce2);names(sce_list)=c(query_species,ifelse(query_species == "Hs", 'Mm', 'Hs'))
      return(sce_list)

    }
    message('--finished--')
  }}

