# 简介
本文档完全基于Seurat V5的相关tutorial。

链接：https://satijalab.org/seurat/articles/get_started_v5_new

时间：2024年11月1日

这一文档主要是运用10 X Genomics提供的2700个PBMCs的数据，一个基础的大概浏览seurat，和一个简单的普通的分析流程的介绍

数据：1. 可以直接使用这个10X的pbmcs的数据 2.建议可以直接使用自己的一些数据进行相关的分析。

数据获取：https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

下载上述文件之后要记得解压缩

# load library
library(dplyr)
library(Seurat)
library(patchwork)

如果出现问题的话，可以问问ai怎么安装

install.packages("Seurat") 一般来说这些比较常见的包是可以直接安装的

# load data
```
## 一般来说上游标准处理的数据，一般是10x，但是近年来也一些文件也变成了h5 文件格式的

pbmc.data <- Read10X(data.dir = "/替换为你文件的下载路径/filtered_gene_bc_matrices/hg19/，（最好和你的R的工作环境处于同一位置）") #如果是10x格式直接这样打开，
getwd() #可以通过这种方式查看你目前R设置的工作环境
setwd() #可以移动你的工作环境

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) #将原来的文件转换为一个seurat格式的文件
pbmc  #简单查看一下这个对象
#看看是不是seurat格式的数据，因为下面的一些函数都要求对象必须是seurat数据
class(pbmc)

```

```
 An object of class Seurat 
 13714 features across 2700 samples within 1 assay 
 Active assay: RNA (13714 features, 0 variable features)
 1 layer present: counts
```


```
Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#运行CreateSeuratObject这个过程的时候如果出现这个问题，无需在意，这是原数据的问题不是我们的问题啊。
```
如果有这些数据的存在就说明这个seurat对象已经成功建立。

接着简单查看一下这个数据的格式（事实上我也是很后面才慢慢理解这个格式的问题的，我觉得暂时不理解也没什么问题。）

![@assay中的结构](https://github.com/Soft283/bioinformatic_learning-progress/blob/Seurat-tutorial/Images/1.png)]

@assay 当中存放的是数据的相关信息

seurat的格式其实有点像是文件夹

@assay$RNA这个文件夹中的存储的RNA的一些表达量的矩阵

我们可以简单的看一下这些数据

```
head(pbmc@assays$RNA@counts)
```

!["dgCMatrix"](https://github.com/Soft283/bioinformatic_learning-progress/blob/Seurat-tutorial/Images/dgCMatrix.png)

虽然不是特别好，但是也足够我们了解这个结构了

可以看到在counts这个矩阵里存储的是基因的表达量，每一行其实是每个基因在不同的细胞当中的表达量。

左边的列名是不同的基因，横着的点因为名字太长了没被显示出来AAACATCA-1,这些其实是不同细胞的标签。

同样类似的，counts这个矩阵一般存放的是原始的基因表达量的数据，data中存放的是normaliztion后的数据，scale.data中存放的是scaledata后的数据。

此外，一般来说seurat的内置函数都是有相应要求的，比如如果你希望对数据进行normaliztion的操作，函数会自动调用counts当中的数据进行相关的计算，包括后续的降维的相关操作，必须使用scaledata的数据。



# standard pre-processing workflow

## QC and selecting cells for further analysis

这一步一般是过滤掉一些表达不太好的数据，QC的指标一般有
1. 每个细胞中基因的数量
  - 低质量细胞和空液滴一班基因数目非常少
  - 可能会存在一些，没有被分开的细胞Cell doublets，中间可能会有异常高的基因计数
2. 细胞内检测到的分子数
  - 线粒体的读段（reads）百分比过高（低质量或者垂死的细胞中线粒体往往表达的指标是比较高的）

一般来说，如果你用的是别人文章里的数据的话，在method的RNAseq的部分都会清楚的描述，需要过滤那些指标
```
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")    #pbmc的自带计算成分一个函数，因为mt基因往往在前面会有一个标志，同样的你也可以用类似的方法来标注，整个测序数据当中的不同来源的细胞，一般测序公司也会给你打上标签
```
```
#查看前几行数据
head(pbmc[["percent.mt"]])
#对应的结果
                 percent.mt
AAACATACAACCAC-1  3.0177759
AAACATTGAGCTAC-1  3.7935958
AAACATTGATCAGC-1  0.8897363
AAACCGTGCTTCCG-1  1.7430845
AAACCGTGTATGCG-1  1.2244898
AAACGCACTGGTAC-1  1.6643551
```

如果你用的是自己的数据的话，可能就需要通过一些可视化来看看，这些过滤指标在整个群体里是怎么样的。
```
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![Vlnplot](https://github.com/Soft283/bioinformatic_learning-progress/blob/Seurat-tutorial/Images/vlnplot.png)

过滤掉上述的数据

```
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

```
#查看数据是否被过滤掉了
pbmc
#对应的结果
An object of class Seurat 
13714 features across 2638 samples within 1 assay 
Active assay: RNA (13714 features, 0 variable features)
```
可以看到前面应该是2700个数据，在这里被过滤掉了一部分。

## Normalizing the data
```
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

对数据进行标准化，这种标准化是最广泛的用于scRNA-seq,让不同的细胞批次具有可比性。

~~哈哈哈我有一个很有意思，关于为什么要标准化的表格，晚点可以给xdm看看~~



## Identification of highly variable features

在这里教程里的话是取出一部分高变基因，之后再去做的下游的降维，但这个顺序还是看大家自己的取舍，毕竟生信归根到底还是大数据分析，怎么选出你需要用的数据得到最好的结果是你一直要思考的内容。
```
#找到高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#定义十个变化最大的基因
top10 <- head(VariableFeatures(pbmc), 10)
#画出点图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

请注意一个问题，这里的高变基因和差异基因是有区别的，其实是完全不一样的部分。

高变基因，是在不同的细胞间里变化差距很大的基因。

差异基因，是在一定条件下（如疾病状态与健康状态、处理组与对照组）表达水平显著变化的基因。


# Scaling the data
缩放数据，为了让每个基因的表达量统一，方便下游的pca分析等过程，不让基因的表达增多影响权重。
```
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

# Perform linear dimensional reduction
对数据进行线性降维，PCA，仅前面提到的高变基因，如果你想用别的数据，请一定要执行scale data这一步
线性降维：在我的理解里，你可以把他直接理解为一个打分，比如你想买一个东西电脑，你可能会从价格，外观，配置等多个方面来衡量
这个电脑到底值不值得买，可能就可以给他打一个分，pca大致就是这样进行的，他会对你的数据差异最明显的几个方面进行打分，把差异越明显的部分叫做pc1（维度1），pc2（维度2）。

```
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```
```
#简单看一下pca降维的前两个维度
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#点图可视化
DimPlot(pbmc, reduction = "pca") + NoLegend()
#热图的可视化
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

```

# Determine the ‘dimensionality’ of the dataset
选取数据的维度，一般来说越到后面的数据维度，他可能数据之间的差异就越不明显，所以完全可以只选择几个数据维度进行下游的分析，那怎么选呢？

可以通过上面的热图，你觉得差异已经不明显的时候就可以不做了

也可以通过下面的
```
ElbowPlot(pbmc)
```
当你看见细胞的变化越来越小变成一条直线的时候就可以不用继续选择了

当然是识别数据的有效的维度，是比较困难的，你可能在接下来的细胞分群当中分的不是特别的好，可能会回到这一步再多加几个数据集

比如在这里我觉得7~12都是可以的，但是也要考虑你的分群的想法，如果你有一定的生物学背景，或者说是对cellmarker很敏感的话，你可能会发现
PC12和13中的MZB 1是一个很经典的 plasmacytoid DCs的标志物，当然具体选择那些维度，需要你的进一步选择。然而，这些群体非常罕见，如果你没有一些知识的储备，很难把他和数据集中的背景噪音分开。所以这一步是非常难的，可能需要你去阅读别人发的论文当中选择的一些marker和你需要的分群的细胞。

在这里，尽量可以稍微多选几个数量pc进行下游分析，（10、15甚至50！）



# Cluster the cells

对细胞进行一个集群,他把一些具有相同特征的细胞放在一起，越接近，越相似。

```
#这里输入的是你前面选择的维度
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#这里选择的是分辨率，分辨率越高，他分的越细，你可能手动命名会越复杂，但是也会相应的让你一个群体变的越来越干净
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

```
#简单看一下分群的情况
head(Idents(pbmc), 5)
```


# Run non-linear dimensional reduction (UMAP/tSNE)
非线性降维，将相似的单元，放在同一个平面上，确定他们之间的相对的空间分布

注意降维都是有代价的，可能都忽略了一些原始数据当中的复杂性，在这里，这些算法，只确保了非常相似基因表达的细胞的位置是在一起的，但是更广泛的细胞相关性等，是没有涉及的，不要依据他们离得很近，就做处一些生物学结论。

```
pbmc <- RunUMAP(pbmc, dims = 1:10)
```
```
#简单可视化
DimPlot(pbmc, reduction = "umap")
```

```
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```
一定要记得把降维后的数据保存下来，因为非线性降维其实是一个随机化的过程，可能你多次降维的结果并不一样，你需要把这次的降维数据保存下来才能保证之后能复现。


# identify cell 

接下来要做的事情是，定义这些群体，他们分别是那些细胞。

虽然我们可能在复现论文的时候文章当中会涉及到一些他是通过那些基因来确定细胞的

但是我们这里还是假装什么都不知道从0开始一步一步的探究这些细胞是那些细胞

## 找出每个cluster的差异基因
```
#这里以细胞群体2为例
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
```
```
#查看这些差异基因
head(cluster2.markers, n = 5)
#结果
             p_val avg_log2FC pct.1 pct.2    p_val_adj
 IL32 2.593535e-91  1.3221171 0.949 0.466 3.556774e-87
 LTB  7.994465e-87  1.3450377 0.981 0.644 1.096361e-82
 CD3D 3.922451e-70  1.0562099 0.922 0.433 5.379250e-66
 IL7R 1.130870e-66  1.4256944 0.748 0.327 1.550876e-62
 LDHB 4.082189e-65  0.9765875 0.953 0.614 5.598314e-61
```
在这里你就得到了这个cluster比较显著的差异基因，这个时候你可以把它丢给ai，或者是用一些网站来确定这些细胞的身份，

```
#一步找出cluster的差异基因，同时只保留那些，变化量大于1的上调基因。
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
```

当然这里面的t检验提供了很多不同的方式，应对不同的计算，在这里不一一演示了

### 差异基因的可视化
```
#小提琴图
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
#FeaturePlot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
#热图
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#点图
DotPlot(pbmc,features = c("MS4A1", "CD79A"))
```

## 通过分子marker确认这些细胞

|Cluster ID	|Markers	|Cell Type|
|-----------|--------|---------|
|0	|IL7R, CCR7 |	Naive CD4+ T|
|1	|CD14, LYZ	| CD14+ Mono|
|2	|IL7R, S100A4 |	Memory CD4+|
|3	|MS4A1	|B|
|4	|CD8A	|CD8+ T|
|5	|FCGR3A, MS4A7	|FCGR3A+ Mono|
|6	|GNLY, NKG7 |	NK|
|7	|FCER1A, CST3	|DC|
|8	|PPBP	|Platelet|

```
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```


```
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
```



