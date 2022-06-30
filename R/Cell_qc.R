#' Cell_qc
#'
#' @param expression_profile
#' @param species
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
Cell_qc<-function(expression_profile,species,metadata) {
  require(Seurat)
  if (dim(expression_profile)[1]==0 ||dim(expression_profile)[2]==0) {
    stop("`expression_profile`should be a cell-gene matrix")
  }
  seob=CreateSeuratObject(counts =expression_profile, min.cells = 3, min.features = 200)
  seob@meta.data=metadata
  if (species=='Mm') {
    seob$percent.mt=PercentageFeatureSet(seob,pattern = '^Mt')
  }
  else if(species=="Hs"){
    seob$percent.mt=PercentageFeatureSet(seob,pattern = '^MT')
  }
  seob<-subset(seob,percent.mt <25)
  seob <- NormalizeData(seob,
                        normalization.method = "LogNormalize")
  if (dim(seob)[2]>=20000) {
    seob<-seob[,sample(dim(seob)[2],20000)]
  }
  return(seob)}



