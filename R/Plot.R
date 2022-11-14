#' Sankyplot
#'
#' @param auc
#' correlation bewtween query and ref celltypes(results from RUN_dscBLAST)
#' @param top_n
#' top n correlated cell_type
#' default: 3
#'
#' @param save
#' default: T
#'
#' @param color
#' vectors of colors
#' default: c("#ffcc99","#66cccc")

#' @param use_shortname
#' default: T
#'
#' @param highlight
#' default: F
#' @param custom.row 
#' default: NULL. If provided, filter rows in the auc matrix.
#' @param custom.col
#' default: NULL. If provided, filter columns in the auc matrix.
#' @return p
#' a sankey plot
#' @export
#'
#' @importFrom reshape2 colsplit
#' @importFrom networkD3 sankeyNetwork saveNetwork
#' @importFrom webshot webshot
#' @importFrom grDevices pdf dev.off
#' @importFrom dplyr %>%
Sanky_plot<- function(auc,top_n=3,color=c("#ffcc99","#66cccc"),use_shortname=T,save=T,cutoff=0.8,highlight=F,custom.row=NULL,custom.col=NULL) {
  if(highlight){auc=auc[['auc_highlight']] %>%as.matrix()
  }else{auc=auc[['auc_total']] %>%as.matrix()}
  rownames(auc)=colsplit(rownames(auc),'_',names = c('c1','c2'))$c2
  env_corr_list <- melt(auc)
  env_corr_list<-env_corr_list[,c(2,1,3)]
  colnames(env_corr_list )<-c("query","ref","AUROC")
  query_name=unique(env_corr_list$query)
  edges <- data.frame(matrix(nrow = 0,ncol = 3))
  colnames(edges)=colnames(env_corr_list)
  for (i in 1:length(query_name)) {
    Row <- which(query_name[i] ==  env_corr_list[,"query", drop = TRUE])
    edges1=env_corr_list[Row,]
    edges1=edges1[sort(edges1$AUROC,decreasing = T,index.return=T)$ix[1:top_n],]
    edges <- rbind(edges,edges1)
  }
  edges=edges[edges$AUROC>=cutoff,]
  
  #filter query and ref celltype
  if(!is.null(custom.row)){edges=edges[edges$ref %in% custom.row,]}
  if(!is.null(custom.col)){edges=edges[edges$query %in% custom.col,]}
  
  d3links <- edges
  d3nodes <- data.frame(name = unique(c(edges$query, edges$ref)), stringsAsFactors = FALSE)
  d3nodes$seq <- 0:(nrow(d3nodes) - 1)
  query_num=dim(d3nodes[grepl('Query',d3nodes$name),])[1]
  ref_num=dim(d3nodes)[1]-query_num
  d3nodes$Family=c(rep('Query',query_num),rep('Ref',ref_num))

  d3links <- merge(d3links, d3nodes, by.x="query", by.y="name")
  names(d3links)[4] <- "source"
  d3links <- merge(d3links, d3nodes, by.x="ref", by.y="name")
  names(d3links)[6] <- "target"
  names(d3links)[3] <- "AUROC"
  d3links <- subset(d3links, select=c("source", "target", "AUROC",'Family.x'))
  colnames(d3links)[4]='Family'
  d3nodes <- subset(d3nodes, select=c("name",'Family'))
  d3nodes$shortname=colsplit(d3nodes$name,'[|]',names = c('c1','c2'))$c2

  cols <- data.frame(
    domain = c("Query","Ref"),
    color = color)
  cols$domain <- sprintf('"%s"', cols$domain)
  cols$color <- sprintf('"%s"', cols$color)


  my_color <- c('d3.scaleOrdinal().domain([',
                paste(c(cols$domain), collapse = ", "),
                "]) .range([",
                paste(c(cols$color), collapse = ", "),
                "])")
  my_color <- paste(my_color, collapse = "")

  if(use_shortname){
    p<-sankeyNetwork(Links = d3links,
                       Nodes = d3nodes,
                       Source = "source",
                       Target = "target",
                       Value = "AUROC",
                       NodeID = "shortname",
                       NodeGroup = "Family",
                       colourScale = my_color,
                       units = "votes",
                     height = length(d3links$target)*30,
                       width =500,
                       nodeWidth = 30,
                       nodePadding = 10,
                       sinksRight = FALSE,
                       fontSize = 12)
    #saving
    if(save){
      message('--Saving Plots--')
      saveNetwork(p,'./sankey.html')
      webshot('./sankey.html','./sankey.pdf')
       system('rm ./sankey.html')
    }

  }
  else{
  p<-sankeyNetwork(Links = d3links,
                     Nodes = d3nodes,
                     Source = "source",
                     Target = "target",
                     Value = "AUROC",
                     NodeID = "name",
                     NodeGroup = "Family",
                     colourScale = my_color,
                     units = "votes",
                     height = length(d3links$target)*30,
                     width =500,
                     nodeWidth = 30,
                     nodePadding = 10,
                     sinksRight = FALSE,
                     fontSize = 12)
  #saving
  if(save){
    message('--Saving Plots--')
    saveNetwork(p,'./sankey.html')
    webshot('./sankey.html','./sankey.pdf')
    system('rm ./sankey.html')
  }
  }
  return(p)

}

#' Networkplot
#' @param auc
#' correlation bewtween query and ref celltypes(results from RUN_dscBLAST)
#' @param top_n
#' top n correlated cell_type
#' default: 3
#' @param save
#' default: F
#'
#' @param color
#' vectors of color
#' default: c("#ffcc99","#66cccc")
#' @param use_shortname
#' default: T
#' @param set.seed
#' random value
#' default: 111
#' @param highlight
#' default: F
#' @param custom.row 
#' default: NULL. If provided, filter rows in the auc matrix.
#' @param custom.col
#' default: NULL. If provided, filter columns in the auc matrix.
#' @export
#' @importFrom reshape2 melt colsplit
#' @importFrom igraph graph_from_data_frame V E
#' @importFrom dplyr %>%
Network_plot<- function(auc,top_n=3,color=c("#ffcc99","#66cccc"),use_shortname=T,save=F,cutoff=0.8,set.seed=111,highlight=F,custom.row=NULL,custom.col=NULL) {
  require(igraph)
  if(highlight){auc=auc[['auc_highlight']] %>%as.matrix()}
  else{auc=auc[['auc_total']] %>%as.matrix()}
  rownames(auc)=colsplit(rownames(auc),'_',names = c('c1','c2'))$c2
  env_corr_list <- melt(auc)
  env_corr_list<-env_corr_list[,c(2,1,3)]
  colnames(env_corr_list )<-c("query","ref","AUROC")
  query_name=unique(env_corr_list$query)
  edges <- data.frame(matrix(nrow = 0,ncol = 3))
  colnames(edges)=colnames(env_corr_list)
  for (i in 1:length(query_name)) {
    Row <- which(query_name[i] ==  env_corr_list[,"query", drop = TRUE])
    edges1=env_corr_list[Row,]
    edges1=edges1[sort(edges1$AUROC,decreasing = T,index.return=T)$ix[1:top_n],]
    edges <- rbind(edges,edges1)
  }
  edges=edges[edges$AUROC>=cutoff,]
  
  #filter query and ref celltype
  if(!is.null(custom.row)){edges=edges[edges$ref %in% custom.row,]}
  if(!is.null(custom.col)){edges=edges[edges$query %in% custom.col,]}
  
  net.igraph = graph_from_data_frame(edges,directed = F)
  V(net.igraph)$color = c(rep(color[1],time=length(unique(edges$query))),rep(color[2],time=length(unique(edges$ref))))
  shortName= colsplit(names(V(net.igraph)),'[|]',names = c('c1','c2'))$c2
  V(net.igraph)$shortName = shortName
  V(net.igraph)$Name =names(V(net.igraph))
  E(net.igraph)$weight = edges$AUROC
  set.seed(set.seed)
  if(use_shortname){
  plot(net.igraph,
       vertex.color = V(net.igraph)$color,
       vertex.frame.color = 'white',
       vertex.label.family = 'sans',
       vertex.label = V(net.igraph)$shortName,
       vertex.label.cex = 0.9,
       vertex.label.color = 'black',
       vertex.label.dist = 0,
       edge.color = 'gray70',
       edge.sizes = E(net.igraph)$weight)
    if(save){
      message('--Saving Plots--')
      pdf('Network.pdf',height = 10,width = 10)
      plot(net.igraph,
           vertex.color = V(net.igraph)$color, 
           vertex.frame.color = 'white',vertex.label.family = 'sans',
           vertex.label = V(net.igraph)$shortName,#use_shortName = T
           vertex.label.cex = 0.9,
           vertex.label.color = 'black',
           vertex.label.dist = 0, 
           edge.sizes = E(net.igraph)$weight)
      dev.off()

    }}
  else{
    plot(net.igraph,
         vertex.color = V(net.igraph)$color,
         vertex.frame.color = 'white',
         vertex.label.family = 'sans',
         vertex.label = V(net.igraph)$Name,
         vertex.label.cex = 0.9,
         vertex.label.color = 'black',
         vertex.label.dist = 0,
         edge.color = 'gray70',
         edge.sizes = E(net.igraph)$weight)

    if(save){
      message('--Saving Plots--')
      pdf('Network.pdf',height = 10,width = 10)
      plot(net.igraph,
           vertex.color = V(net.igraph)$color,
           vertex.frame.color = 'white',
           vertex.label = V(net.igraph)$shortName,
           vertex.label.cex = 0.9,
           vertex.label.color = 'black',
           vertex.label.dist = 0, 
           edge.sizes = E(net.igraph)$weight)
      dev.off()

    }}

}

#' Heatmap plot
#'
#' @param auc
#' correlation bewtween query and ref celltypes(results from RUN_dscBLAST)
#' @param top_n
#' top n correlated cell_type
#' default: 3
#'
#' @param save
#' default: F
#'
#' @param color
#' vectors of color
#' default:colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
#' @param use_shortname
#' default: T
#' @param refine
#' default: T
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape2 colsplit melt
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr %>%
Heatmap_plot<- function(auc,top_n=3,use_shortname=T,color=colorRampPalette(brewer.pal(n = 7,name = "GnBu"))(100),save=F,cutoff=0.8,refine=T,highlight=F,custom.row=NULL,custom.col=NULL) 
  {
  if (highlight) {
    auc = auc[["auc_highlight"]] %>% as.matrix()
  }
  else {
    auc = auc[["auc_total"]] %>% as.matrix()
  }
  # rownames(auc) = colsplit(rownames(auc), "_", names = c("c1", 
  #                                                        "c2"))$c2
  env_corr_list <- melt(auc)
  env_corr_list <- env_corr_list[, c(2, 1, 3)]
  colnames(env_corr_list) <- c("query", "ref", "AUROC")
  query_name = unique(env_corr_list$query)
  edges <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(edges) = colnames(env_corr_list)
  for (i in 1:length(query_name)) {
    Row <- which(query_name[i] == env_corr_list[, "query", 
                                                drop = TRUE])
    edges1 = env_corr_list[Row, ]
    edges1 = edges1[sort(edges1$AUROC, decreasing = T, index.return = T)$ix[1:top_n], 
    ]
    edges <- rbind(edges, edges1)
  }
  edges = edges[edges$AUROC >= cutoff, ]
  edges$ref=colsplit(edges$ref,'_',names = c('c1','c2'))$c2
  if (!is.null(custom.row)) {
    edges = edges[edges$ref %in% custom.row, ]
  }
  if (!is.null(custom.col)) {
    edges = edges[edges$query %in% custom.col, ]
  }
  select = unique(edges$ref) %>% as.character()
  select_query = unique(edges$query) %>% as.character()
  use_heatmap = auc[select, select_query] %>% as.data.frame()
  if (dim(use_heatmap)[2] > 1) {
    if (refine & is.null(custom.row) & is.null(custom.col)) {
      for (i in 1:dim(use_heatmap)[2]) {
        use_heatmap = use_heatmap[order(-use_heatmap[[i]]), 
        ]
        use_heatmap[, i][(top_n + 1):dim(use_heatmap)[1]] = 0
      }
    }
    use_heatmap = use_heatmap[select, select_query] %>% as.matrix()
    colnames(use_heatmap) = select_query
  }else {
    colnames(use_heatmap) = select_query
    rownames(use_heatmap) = select
  }
  use_heatmap[use_heatmap < cutoff] = 0
  if (use_shortname) {
    rownames(use_heatmap) = colsplit(rownames(use_heatmap), 
                                     "[|]", names = c("c1", "c2"))$c2
    pheatmap(use_heatmap, border = "white", legend = T, legend_breaks = c(0, 
                                                                          0.5, 1), treeheight_row = 20, treeheight_col = 20, 
             cluster_cols = F, cluster_rows = F, color = color)
    if (save) {
      message("--Saving Plots--")
      pdf("./heatmap.pdf", height = 10, width = 10)
      pheatmap(use_heatmap, border = "white", legend = T, 
               legend_breaks = c(0, 0.5, 1), treeheight_row = 20, 
               treeheight_col = 20, cluster_cols = F, cluster_rows = F, 
               color = color)
      dev.off()
    }
  }
  else {
    pheatmap(use_heatmap, border = "white", legend = T, legend_breaks = c(0, 
                                                                          0.5, 1), treeheight_row = 20, treeheight_col = 20, 
             cluster_cols = F, cluster_rows = F, color = color)
    if (save) {
      message("--Saving Plots--")
      pdf("./heatmap.pdf", height = 10, width = 10)
      pheatmap(use_heatmap, border = "white", legend = T, 
               legend_breaks = c(0, 0.5, 1), treeheight_row = 20, 
               treeheight_col = 20, cluster_cols = F, cluster_rows = F, 
               color = color)
      dev.off()
    }
  }
}

#' Plot markers
#'
#' @param sce
#' a single cell experiment object
#' @param species
#' options: c('Hs','Mm')
#'
#' @param features
#' genes offered for visualization
#' @export
#' @importFrom reshape2 melt
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot 
plotMarkers <- function(sce,species,features,color=NULL){
  require(ggplot2)
  if(is.list(sce)){sce=sce[[species]]}
  gene=rownames(sce@assays@data$counts)
  data=sce@assays@data$counts
  "%ni%" <- Negate("%in%")
  for (i in features) {
    if (i %ni% gene){
      stop(i," not found ")}
  }
  vln.df=as.data.frame(data[features,])
  vln.df$gene=rownames(vln.df)
  vln.df=melt(vln.df,id="gene")
  colnames(vln.df)[c(2,3)]=c("CB","exp")

  anno=data.frame(CB=colnames(sce),celltype=sce$cell_type)
  vln.df=inner_join(vln.df,anno,by="CB")
  #color choose
  mycolor = c(c("#db6968","#4d97cd","#99cbeb","#459943",
                "#fdc58f","#e8c559","#a3d393","#f8984e"),
              brewer.pal(11,'Spectral'),
              brewer.pal(11,'PRGn')[c(1:4)],
              brewer.pal(11,'BrBG')[c(7:11,1:5)],
              brewer.pal(11,'PiYG')[c(4:1)],
              brewer.pal(11,'RdBu')[c(3:5,7:11)],
              brewer.pal(11,'PiYG')[c(11:7)],
              brewer.pal(11,'RdGy')[c(7:10)],
              brewer.pal(12,'Set3')[c(3:12)]
  )
  if(is.null(color)){color=colorRampPalette(mycolor)(length(features))}
  #plot
  p=vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=gene),scale = "width")+
    facet_grid(vln.df$gene~.,scales = "free_y")+
    scale_fill_manual(values = color)+
    scale_x_discrete("")+ylab('Expression level')+
    theme_bw()+
    theme(
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  return(p)
}
