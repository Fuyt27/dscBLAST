#' Gene transfer
#'
#' @param expression_profile
#' cell-gene matrix
#' @param species
#'
#' @return
#' @export
#'
#' @examples
transfer_gene <- function(expression_profile,species){
  if(species=='Hs'){
    Mouse_gene <- gene_transfer[gene_transfer$Human_gene%in%rownames(expression_profile),]
    Mouse_gene<-Mouse_gene[!duplicated(Mouse_gene$Mouse_gene),]
    Mouse_gene<-Mouse_gene[!Mouse_gene$Mouse_gene=='',]
    expression_profile <- expression_profile[Mouse_gene$Human_gene,]
    rownames(expression_profile) <-Mouse_gene$Mouse_gene
    return(expression_profile)
  }
  if(species=='Mm'){
    Human_gene <- gene_transfer[gene_transfer$Mouse_gene%in%rownames(expression_profile),]
    Human_gene<-Human_gene[!duplicated(Human_gene$Human_gene),]
    Human_gene<-Human_gene[!Human_gene$Human_gene=='',]
    expression_profile <- expression_profile[Human_gene$Mouse_gene,]
    rownames(expression_profile) <-Human_gene$Human_gene
    return(expression_profile)
  }
}
