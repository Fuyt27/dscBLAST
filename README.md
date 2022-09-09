# README

dscBLAST: developing single cell Basic Local Alignment Search Tool

We collected more than 160 datasets (Organ_First author) form approximately 50 articles. These data were grouped into five references: Ectoderm, Endoderm, Mesoderm, pre-organogenesis and extra-embryonic tissues. Then MetaNeighbor was used to calculate the similarity between query cells and reference cells. The combination of comprehensive references and MetaNeighbor enables the users to confirm the cell identity and speculate its possible developmental stages rapidly and accurately. dscBLAST is particularly better at developing data, and it also performs very well on adult data.


## Workflow
![](Figure/workflow.jpg)

## **Download**

```shell
wget https://github.com/Fuyt27/dscBLAST/archive/refs/tags/1.0.1.tar.gz
```

## **Install**
```R
install.package('dscBLAST_1.0.1.tar.gz',repos=NULL,type='source')
```

## **Tutorial**

For demonstration purpose, we took the dataset **[ectoderm](https://github.com/Fuyt27/dscBLAST/releases/download/1.0.1/ectoderm.rds)** as example.

#### Package loading

```R
#make sure the following packages are installed in your R-library.
library(MetaNeighbor)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(ggplot2)
library(igraph)
library(networkD3)
library(webshot)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(stringr)
```

#### Data processing

- detect MT genes
- down-sample
  - default: 20000 cells
- normalization
  - default: SCT transform
  - if a normalized matrix is offered, set  `mtx.type`='normalized' [default: 'raw']
- create dscBLAST object

```R
library(dscBLAST)
#input data loading
seob<-readRDS('./ectoderm.rds')
count<-seob@assays$RNA@counts #input gene expression matrix
meta<-seob@meta.data #input metadata(set columname:'cell_type') 

#data processing
##generally we run dscblast in both human and mouse datasets(recommanded!), if you prefer to run dscblast for one species, please set ref_species ='single'.
## if you want to normalize your data depending on batch, please set `batch` ='your_batchname'.
## if the expression_profile offered is normalized in advance, please set `mtx.type` = 'normalized'.
sce<-create_dscBLASTobject(expression_profile = count,query_species= 'Hs',metadata = meta,downsample = 20000,batch ='default',mtx.type = 'raw',ref_species ='both')
```



#### Run dscBLAST

```R
#run dscBLAST
auc <- RUN_dscBLAST(sce,query_species = 'Hs',reference = 'all',highlight_celltype=NULL)

##generally we run dscblast in all datasets, if you prefer to run dscblast using certain datasets or datasets, please change the option `reference` ,for example:
auc <- RUN_dscBLAST(sce,query_species = 'Hs',reference = c('ecto','pre'),highlight_celltype=NULL)

## if you want to highlight certain cell types, please set `highlight_celltype`  
auc_highlight <- RUN_dscBLAST(sce,query_species = 'Hs',reference = c('ecto','pre'),highlight_celltype='Melanocyte')
```

#### Check auc
```R
#get total auc 
auc_result <- auc[['auc_total']]

#get highlight auc 
auc_highlight_result <- auc_highlight[['auc_highlight']]
```


#### Visualization

we provide three forms of visualization

```R
#Sankey plot
Sanky_plot(auc)

#Network plot
Network_plot(auc)

#Heatmap plot
Heatmap_plot(auc)

#if you want to show highlight cell types only,please set `highlight=T`
Sanky_plot(auc_highlight,highlight=T)
```



#### Marker

```R
# help to get markers from certain datasets and cell types,for example:
marker_list <- getMarkers(species ='Hs',reference ='ecto' ,ref_celltype ='Melanocyte',ref_dataset ='Eye_Han' ,marker_top_n =20)

#Plot markers in your data
plotMarkers(sce,species = 'Hs',features = c('MLANA','TYRP1')) #Melanocyte markers
```



#### Run dscBLAST with stage info
make sure you have downloaded  **[ref.zip](https://github.com/Fuyt27/dscBLAST/releases/download/1.0.1/ref.zip)** and save it in your own local directory
```R
#run dscBLAST with stage info
auc2 <- RUN_dscBLAST_stage(sce,auc,query_species = 'Hs',ref_dir = 'your_local_dir')

#visualization
Heatmap_plot(auc2)
Network_plot(auc2)
Sanky_plot(auc2)
```



## Citation



## Contact

daiyuanliu@zju.edu.cn
