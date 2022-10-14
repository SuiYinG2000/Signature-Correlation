# 基因Signature与GO本体论和通路之间的相关性分析  
## 简介  
此代码是较为用户友好型的脚本，用于**某一特征基因集**与**生物过程、通路**之间的**相关性分析**，其中可供选择的参数有很多，可以自定义参数组合进行相关性分析。  
PS：也可以用于**单个基因**与生物过程、通路之间的相关性分析。  
此外还提供了将GeneID的基因集统一转换为Symbol以及批量运行多种参数组合的批处理脚本。  
## 下载与安装  
直接下载并在本地使用  
## 内容  
    * GeneSignature-Pathways_Correlation.Rmd：相关性分析  
    * GSEA_GeneSets_ID2Symbol.Rmd：GeneID基因集转换为Symbol  
    * GSEA_GeneSets_Normalize2Symbol.Rmd：GeneName基因集转换为Symbol的脚本（不推荐
    * loop.R：多种参数组合的批量相关性分析
    * .\lib\GeneSignature-Pathways_Correlation_Function.R：GeneSignature-Pathways_Correlation.Rmd的配套自定义函数  
## 使用  
### GeneSignature-Pathways_Correlation.Rmd  

1、 材料  
注意：运行脚本需要将所有原始数据放置于工作目录下的data文件夹中  
    1.1、 表达谱数据：来自TCGA或CGGA的数据，或者其他表达谱数据  
        数据要求：行为基因，列为样本，且为csv格式  
        1）TCGA数据：可以从TCGA官网下载数据并通过script：TCGA_Rawdata_Merge-Expression.Rmd，或者其他方式获得  
        2）CGGA数据：可以从CGGA官网下载数据并通过Project：[CGGA_RawData_Merge-Combat](https://github.com/SuiYinG2000/CGGA_RawData_Merge-Combat)，或其他方式获得  
        3）其他表达谱数据  
    1.2、 Gene名称对应Symbol注释文件  
        推荐从NCBI GENE上下载[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    1.3、 GO本体论和通路基因集  
        数据要求：需为gmt格式  
        推荐从GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)网站上下载**Gene ID**组成的基因集文件  
        （原因：若下载Gene Symbols的文件，基因集会少于网站所展示的数量，同时部分基因并不是Symbol）  
        需要同时下载**GO本体论**和**Pathway**数据  
        **注意**可以通过**GSEA_GeneSets_ID2Symbol.Rmd**或者**GSEA_GeneSets_Normalize2Symbol.Rmd**将基因集统一转换为Symbol  
    1.4、配套自定义函数集  
        .\lib\GeneSignature-Pathways_Correlation_Function.R：GeneSignature-Pathways_Correlation.Rmd  

2、 参数  
    * `workdir`：工作目录  
    * `libdir`：配套自定义函数集  
    * `expression_file`：表达谱数据  
        >>* "\\data\\XXX.csv"    
    * `GOGeneSet`：GO基因集  
        * "\\data\\XXX.gmt"    
    * `PWGeneSet`：Pathway基因集  
        * "\\data\\XXX.gmt"    
    * `GeneSet_Name_Normalized`：基因集基因名称是否通过GSEA_GeneSets_ID2Symbol.Rmd或GSEA_GeneSets_Normalize2Symbol.Rmd进行转换  
        * TRUE  
        * FALSE    
    * `GeneNameAnno`：Gene名称对应Symbol注释文件  
        * "\\data\\XXX.csv"    
    * `SignatureGenes_Name`：基因Signature的名称或单个基因名称  
        * "m6A_Signature_Genes"  
        * "STK26"  
    * `SignatureGenes`：基因Signature中基因集，用逗号","分隔  
    * `corAnalysis_Method`：相关性分析方法：  
        * "pearson"  
        * "spearman"  
    * `calOrder`：计算Signature时标准化的选择：  
        * standFirst：每个基因先在所有样本中进行标准化，再计算Signature Score  
        * calFirst：先计算Signature Score，再进行标准化  
        * noStand：不进行标准化  
    * `pickMethod`：计算Signature Score的方法：  
        * mean：每个样本中Signature Genes的平均表达值  
        * median：每个样本中Signature Genes的中间表达值  
    * `log`：对数转换  
        * log2(n)  
        * log2(n+1)  
        * FALSE  
    * `exp`：取幂  
        * TRUE  
        * FALSE  


