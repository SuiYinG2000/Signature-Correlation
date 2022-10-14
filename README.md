# 基因Signature与GO本体论和通路之间的相关性分析  
## 简介
此代码是较为用户友好型的脚本，用于分析某一特征基因集与生物过程、通路之间的相关性，其中可供选择的参数有很多，可以自定义参数组合进行相关性分析。  
此外还提供了将GeneID的基因集统一转换为Symbol以及批量运行多种参数组合的批处理脚本。  
## 下载与安装
可以直接下载并在本地使用  
## 内容
    * GeneSignature-Pathways_Correlation.Rmd：相关性分析  
    * GSEA_GeneSets_ID2Symbol.Rmd：GeneID基因集转换为Symbol  
    * GSEA_GeneSets_Normalize2Symbol.Rmd：GeneName基因集转换为Symbol的脚本（不推荐
    * loop.R：多种参数组合的批量相关性分析
    * .\lib\GeneSignature-Pathways_Correlation_Function.R：GeneSignature-Pathways_Correlation.Rmd的配套自定义函数