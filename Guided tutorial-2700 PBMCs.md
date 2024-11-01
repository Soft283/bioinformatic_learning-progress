# 简介
本文档完全基于Seurat V5的相关tutorial。
链接：https://satijalab.org/seurat/articles/get_started_v5_new
时间：2024年11月1日
这一文档主要是运用10 X Genomics提供的2700个PBMCs的数据，一个基础的大概浏览seurat，和一个简单的普通的分析流程的介绍
数据：1. 可以直接使用这个10X的pbmcs的数据 2.建议可以直接使用自己的一些数据进行相关的分析。
数据获取：https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# load data
```
## 一般来说上游标准处理的数据，一般是10x，但是近年来也一些文件也变成了h5 文件格式的
pbmc.data <- Read10X(data.dir = "/替换为你文件的下载路径，（最好和你的R的工作环境处于同一位置）") #如果是10x格式直接这样打开，
getwd() #可以通过这种方式查看你目前R设置的工作环境

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) #将原来的文件转换为一个seurat格式的文件

```

