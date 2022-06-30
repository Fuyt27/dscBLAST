#' @title
#' create dscBLAST object
#'
#' @param expression_profile
#' cell-gene matrix
#' @param species
#' Hs/Mm
#' @param trans
#' boolean value
#' @param metadata
#' meta-data
#' @return Returns a Sce object or list
#'
#'
#' @export
#'
#' @examples
#' sce<-createdscBLASTobject()
#' @importFrom Seurat PercentageFeatureSet as.SingleCellExperiment CreateSeuratObject SCTransform
#' @importFrom SingleCellExperiment
#' @importFrom dplyr inner_join join
create_dscBLASTobject <- function(expression_profile,species,trans,metadata) {
  if(!trans){
    seob=Cell_qc(expression_profile,species,metadata)
    seob<- SCTransform(seob,
                       variable.features.n = 3000,
                       vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
                       verbose = T)
    seob_sct = CreateSeuratObject(counts = seob@assays[["SCT"]]@data)
    seob_sct@meta.data=seob@meta.data
    sce=as.SingleCellExperiment(seob_sct)
    sce$study_id='Query'
    return(sce)}
  if(trans){
    expression_profile_transfer=transfer_gene(expression_profile,species)
    seob=Cell_qc(expression_profile,species,metadata)
    seob<- SCTransform(
      seob,
      variable.features.n = 3000,
      vars.to.regress = c("nCount_RNA", "nFeature_RNA",'percent.mt'),
      verbose = T)

    counts = seob@assays[["SCT"]]@data
    seob_sct = CreateSeuratObject(counts = counts)
    seob_sct@meta.data=seob@meta.data
    sce1=as.SingleCellExperiment(seob_sct)
    sce1$study_id='Query'

    counts2<-transfer_gene(counts,species)
    seob_sct = CreateSeuratObject(counts = counts2)
    seob_sct@meta.data=seob@meta.data
    sce2=as.SingleCellExperiment(seob_sct)
    sce2$study_id='Query'
    return(list(sce1,sce2))
  }}

