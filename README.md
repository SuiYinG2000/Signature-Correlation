# 基因Signature与GO功能注释和通路之间的相关性分析  

## 简介  

此代码是较为用户友好型的脚本，用于**某一特征基因集**与**功能、通路**之间的**相关性分析**，其中有较多可供选择的参数，且可以自定义参数组合进行相关性分析。  
PS：也可以用于**单个基因**与功能、通路之间的相关性分析。  
此外还提供了将基因集统一转换为Symbol以及批量运行多种参数组合的批处理脚本。  

## 下载与安装  

直接下载并在本地使用  

## 内容  

+ GeneSignature-Pathways_Correlation.Rmd：相关性分析  
+ GSEA_GeneSets_ID2Symbol.Rmd：GeneID基因集转换为Symbol  
+ GSEA_GeneSets_Normalize2Symbol.Rmd：GeneName基因集转换为Symbol的脚本（不推荐
+ loop.R：多种参数组合的批量相关性分析
+ .\lib\GeneSignature-Pathways_Correlation_Function.R：GeneSignature-Pathways_Correlation.Rmd的配套自定义函数  

## 介绍与使用  

### 一、GeneSignature-Pathways_Correlation.Rmd  

#### （一）介绍

该代码通过使用基因表达谱信息，分析某一个Gene Signature与GO功能注释、Pathway通路之间的相关性分析。Gene Signature是一组具有生物学意义的基因集合，通过该组基因集的表达计算该Signature的打分，来代表该Signature在样本中的表达情况。  
（ps：支持输入单个基因，计算该基因与各功能、通路之间的相关性

#### （二）使用

1. 材料  
**注意**：运行脚本需要将所有原始数据放置于工作目录下的data文件夹中  
    1. 表达谱数据：来自TCGA或CGGA的数据，或者其他表达谱数据  
    **数据要求**：行为基因，列为样本，且为csv格式  
        + TCGA数据：可以从TCGA官网下载数据并通过script：TCGA_Rawdata_Merge-Expression.Rmd，或者其他方式获得  
        + CGGA数据：可以从CGGA官网下载数据并通过Project：[CGGA_RawData_Merge-Combat](https://github.com/SuiYinG2000/CGGA_RawData_Merge-Combat)，或其他方式获得  
        + 其他表达谱数据  
    2. Gene名称对应Symbol注释文件  
    推荐从NCBI GENE上下载[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    3. GO功能注释和通路基因集  
    **数据要求**：需为gmt格式  
    推荐从GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)网站上下载**Gene ID**组成的基因集文件  
    （原因：若下载Gene Symbols的文件，基因集会少于网站所展示的数量，同时部分基因并不是Symbol）  
    **注意**：
        + 需要同时下载**GO功能注释**和**Pathway**的基因集，或者内部手动修改代码从而仅计算与其中一个的相关性。  
        + 可以通过**GSEA_GeneSets_ID2Symbol.Rmd**或者**GSEA_GeneSets_Normalize2Symbol.Rmd**将基因集统一转换为Symbol  
    4. 配套自定义函数集  
.\lib\GeneSignature-Pathways_Correlation_Function.R：GeneSignature-Pathways_Correlation.Rmd  

2. 参数  
**建议**所有路径信息都为绝对路径。若将该脚本放置于工作目录下双击打开，可以选择使用相对路径
    + `workdir`：工作目录  
    + `libdir`：配套自定义函数集  
    + `expression_file`：表达谱数据  
      + "\\data\\XXX.csv"  
    + `GOGeneSet`：GO基因集
      + "\\data\\XXX.gmt"
    + `PWGeneSet`：Pathway基因集  
      + "\\data\\XXX.gmt"
    + `GeneSet_Name_Normalized`：基因集基因名称是否通过GSEA_GeneSets_ID2Symbol.Rmd或GSEA_GeneSets_Normalize2Symbol.Rmd进行转换  
      + TRUE  
      + FALSE
    + `GeneNameAnno`：Gene名称对应Symbol注释文件  
      + "\\data\\XXX.csv"
    + `SignatureGenes_Name`：基因Signature的名称或单个基因名称  
      + "m6A_Signature_Genes"  
      + "STK26"  
    + `SignatureGenes`：基因Signature中基因集，用逗号","分隔  
    + `corAnalysis_Method`：相关性分析方法：  
      + "pearson"  
      + "spearman"  
    + `calOrder`：计算Signature时标准化的选择：  
      + "standFirst"：每个基因先在所有样本中进行标准化，再计算Signature Score  
      + "calFirst"：先计算Signature Score，再进行标准化  
      + "noStand"：不进行标准化  
    + `pickMethod`：计算Signature Score的方法：  
      + "mean"：每个样本中Signature Genes的平均表达值  
      + "median"：每个样本中Signature Genes的中间表达值  
    + `log`：对数转换  
      + "log2(n)"  
      + "log2(n+1)"  
      + FALSE  
    + `exp`：取幂  
      + TRUE  
      + FALSE  

### 二、GSEA_GeneSets_ID2Symbol.Rmd

#### （一）介绍

该代码用于将从GESA MSigDB网站上下载的**Gene IDs**组成的基因集，统一转换成**Gene Symbol**的形式  

#### （二）使用

1. 材料  
**注意**：运行脚本需要将所有原始数据放置于工作目录下的data文件夹中  
    1. Gene名称对应Symbol注释文件  
    推荐从NCBI GENE上下载[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    2. GO功能注释和通路的GeneID基因集  
    **数据要求**：需为gmt格式  
    从GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)网站上下载**Gene ID**组成的基因集文件  
    **注意**：
        + 需要同时下载**GO功能注释**和**Pathway**的基因集，或者内部手动修改代码从而仅计算与其中一个的相关性。  
2. 参数  
    + `workdir`：工作目录  
    + `libdir`：配套自定义函数集  
    + `resultsdir`：结果文件夹（不建议更改
    + `picsdir`：图片文件夹（不建议更改
    + `GOGeneset`：GO基因集
      + "\\data\\XXX.gmt"
    + `PWGeneSet`：Pathway基因集  
      + "\\data\\XXX.gmt"
    + `GeneNameAnno`：Gene名称对应Symbol注释文件  
      + "\\data\\XXX.csv"

### 三、GSEA_GeneSets_Normalize2Symbol.Rmd

#### （一）介绍

该代码用于将从GESA MSigDB网站上下载的**Gene Symbols**组成的基因集，统一转换成**Gene Symbol**的形式  

#### （二）使用

1. 材料  
**注意**：运行脚本需要将所有原始数据放置于工作目录下的data文件夹中  
    1. Gene名称对应Symbol注释文件  
    推荐从NCBI GENE上下载[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    2. GO功能注释和通路的GeneID基因集  
    **数据要求**：需为gmt格式  
    从GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)网站上下载**Gene Symbols**组成的基因集文件  
    **注意**：
        + 需要同时下载**GO功能注释**和**Pathway**的基因集，或者内部手动修改代码从而仅计算与其中一个的相关性。  
2. 参数  
    + `workdir`：工作目录  
    + `libdir`：配套自定义函数集  
    + `resultsdir`：结果文件夹（不建议更改
    + `picsdir`：图片文件夹（不建议更改
    + `GOGeneset`：GO基因集
      + "\\data\\XXX.gmt"
    + `PWGeneSet`：Pathway基因集  
      + "\\data\\XXX.gmt"
    + `GeneNameAnno`：Gene名称对应Symbol注释文件  
      + "\\data\\XXX.csv"

### 四、loop.R

#### （一）介绍

该脚本为相关性分析的批处理脚本，可以输入不同参数组合进行批量运算。

#### （二）使用

1. 参数
    + `rmd`：GeneSignature-Pathways_Correlation.Rmd所在位置和该rmd文件
    + `workdir`：工作目录
    + `libdir`：配套自定义函数集目录
    + `expr_file`：表达谱数据
      + "\\data\\XXX.csv"
    + `GOGeneset`：GO基因集
      + "\\data\\XXX.gmt"
    + `PWGeneSet`：Pathway基因集  
      + "\\data\\XXX.gmt"
    + `GeneSet_Name_Normalized`：基因集基因名称是否通过GSEA_GeneSets_ID2Symbol.Rmd或GSEA_GeneSets_Normalize2Symbol.Rmd进行转换  
      + TRUE  
      + FALSE
    + `GeneNameAnno`：Gene名称对应Symbol注释文件  
      + "\\data\\XXX.csv"
    + `SignatureGenes_Name`：基因Signature的名称或单个基因名称  
      + "m6A_Signature_Genes"  
      + "STK26"  
    + `SignatureGenes`：基因Signature中基因集，用逗号","分隔  

    + `corAnalysis_Methods`：相关性分析方法（单个字符串/向量）
      + "pearson"  
      + "spearman"  
    + `calOrders`：计算Signature时标准化的选择（单个字符串/向量）
      + "standFirst"：每个基因先在所有样本中进行标准化，再计算Signature Score  
      + "calFirst"：先计算Signature Score，再进行标准化  
      + "noStand"：不进行标准化  
    + `pickMethods`：计算Signature Score的方法（单个字符串/向量）
      + "mean"：每个样本中Signature Genes的平均表达值  
      + "median"：每个样本中Signature Genes的中间表达值  
    + `exps`：取幂（单个字符串/向量）
      + TRUE  
      + FALSE  
    + `logs`：对数转换（单个字符串/向量）
      + "log2(n)"  
      + "log2(n+1)"  
      + FALSE  

## 结语

这个Project是我在硕士生涯中第一份工作，由于参考文献中没有完善的流程介绍，所以一切的一切都是从零开始摸索，包括（TCGA和）CGGA的数据下载及预处理（去批次效应）；包括怎么计算Signature的打分（median/exp的组合是另一篇老师推给的文献中使用的方法；而log(n+1)是GEPIA2中使用的方法；预先进行标准化，并且取平均值是肖肖老师给出的很棒的指导，对此我深信不疑）；包括最后相关性分析的使用。耗时两个月经过艰辛的过程（hh好像也没有特别的艰辛），最终摸索出来的这么一套Protocal。  

截至我完善这篇README文档为止，文献的图也没有能够完全的复现出来的。喜人的是，使用我认为较为正确的一套参数组合后，经TCGA和CGGA数据计算出的强正相关性的功能与通路基本一致，但强负相关性的overlap结果差强人意，CGGA的最终计算得到的负相关性的数据也很少。后续可能要继续尝试其他参数组合（CGGA的样本量太大每次计算一个参数组合就要3、4h左右orz）。  

总之，能够近乎完成这样的一个Project内心多少还是有一点自豪感的，都说代码人最快乐的时候就是当自己的代码普适性很高，于是乎就可以拿着一杯果茶看着代码在后台自己慢慢跑，就好像看着自己长大了的小孩子一样。  

如果你有兴趣看我这比较水的代码并提出一些改良的意见，我会非常感谢你的:)。  

就这样，Peace out.
