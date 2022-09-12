#' Title:RUN_dscBLAST_stage
#'
#' @param sce
#' singlecellexperiment object
#' @param auc
#' results from RUN_dscBLAST
#' @param query_species
#' Hs/Mm
#'
#' @param top_n
#' default: 3
#' @param marker_top_n
#' default: 20
#' @param ref_dir
#' ref data and meta stored here
#'
#' @return auc_final
#' @export
#'
RUN_dscBLAST_stage<- function(sce,auc,query_species,top_n=3,marker_top_n=20,ref_dir,highlight=F) {
  if(is.list(auc) & 'auc_total' %in% names(auc)){auc=auc[['auc_total']]}
  message('--step1:subset data--')
  ref_list=subset_data(auc,top_n,ref_dir)

  message('--step2:subset markers--')
  ref_marker_list=subset_marker(ref_list,marker_top_n)

  message('--step3:train model with stage info--')
  pretrained_model=list()
  for (i in names(ref_list)) {
    var_genes = ref_marker_list[[i]]
    ref_list[[i]] <-ref_list[[i]][var_genes,]
  }
  for (i in names(ref_list)) {
    pretrained_model[[i]] <- trainModel_stage(var_genes = rownames(ref_list[[i]]),dat = ref_list[[i]],study_id = ref_list[[i]]$dataset,cell_type =ref_list[[i]]$cell_type,stage_id=ref_list[[i]]$stage)
    }

  message('--step4:calculate auc using Metaneighbor--')
  auc_final=auc_calculate_withstage(sce,pretrained_model,query_species,highlight)
  message('--finished--')
  return(auc_final)
}
