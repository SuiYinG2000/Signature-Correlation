# ����Signature��GO�����ۺ�ͨ·֮�������Է���  
## ���  
�˴����ǽ�Ϊ�û��Ѻ��͵Ľű�������**ĳһ��������**��**������̡�ͨ·**֮���**����Է���**�����пɹ�ѡ��Ĳ����кܶ࣬�����Զ��������Ͻ�������Է�����  
PS��Ҳ��������**��������**��������̡�ͨ·֮�������Է�����  
���⻹�ṩ�˽�GeneID�Ļ���ͳһת��ΪSymbol�Լ��������ж��ֲ�����ϵ�������ű���  
## �����밲װ  
ֱ�����ز��ڱ���ʹ��  
## ����  
    * GeneSignature-Pathways_Correlation.Rmd������Է���  
    * GSEA_GeneSets_ID2Symbol.Rmd��GeneID����ת��ΪSymbol  
    * GSEA_GeneSets_Normalize2Symbol.Rmd��GeneName����ת��ΪSymbol�Ľű������Ƽ�
    * loop.R�����ֲ�����ϵ���������Է���
    * .\lib\GeneSignature-Pathways_Correlation_Function.R��GeneSignature-Pathways_Correlation.Rmd�������Զ��庯��  
## ʹ��  
### GeneSignature-Pathways_Correlation.Rmd  

1�� ����  
ע�⣺���нű���Ҫ������ԭʼ���ݷ����ڹ���Ŀ¼�µ�data�ļ�����  
    1.1�� ��������ݣ�����TCGA��CGGA�����ݣ������������������  
        ����Ҫ����Ϊ������Ϊ��������Ϊcsv��ʽ  
        1��TCGA���ݣ����Դ�TCGA�����������ݲ�ͨ��script��TCGA_Rawdata_Merge-Expression.Rmd������������ʽ���  
        2��CGGA���ݣ����Դ�CGGA�����������ݲ�ͨ��Project��[CGGA_RawData_Merge-Combat](https://github.com/SuiYinG2000/CGGA_RawData_Merge-Combat)����������ʽ���  
        3���������������  
    1.2�� Gene���ƶ�ӦSymbolע���ļ�  
        �Ƽ���NCBI GENE������[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    1.3�� GO�����ۺ�ͨ·����  
        ����Ҫ����Ϊgmt��ʽ  
        �Ƽ���GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)��վ������**Gene ID**��ɵĻ����ļ�  
        ��ԭ��������Gene Symbols���ļ������򼯻�������վ��չʾ��������ͬʱ���ֻ��򲢲���Symbol��  
        ��Ҫͬʱ����**GO������**��**Pathway**����  
        **ע��**����ͨ��**GSEA_GeneSets_ID2Symbol.Rmd**����**GSEA_GeneSets_Normalize2Symbol.Rmd**������ͳһת��ΪSymbol  
    1.4�������Զ��庯����  
        .\lib\GeneSignature-Pathways_Correlation_Function.R��GeneSignature-Pathways_Correlation.Rmd  

2�� ����  
    * `workdir`������Ŀ¼  
    * `libdir`�������Զ��庯����  
    * `expression_file`�����������  
        >>* "\\data\\XXX.csv"    
    * `GOGeneSet`��GO����  
        * "\\data\\XXX.gmt"    
    * `PWGeneSet`��Pathway����  
        * "\\data\\XXX.gmt"    
    * `GeneSet_Name_Normalized`�����򼯻��������Ƿ�ͨ��GSEA_GeneSets_ID2Symbol.Rmd��GSEA_GeneSets_Normalize2Symbol.Rmd����ת��  
        * TRUE  
        * FALSE    
    * `GeneNameAnno`��Gene���ƶ�ӦSymbolע���ļ�  
        * "\\data\\XXX.csv"    
    * `SignatureGenes_Name`������Signature�����ƻ򵥸���������  
        * "m6A_Signature_Genes"  
        * "STK26"  
    * `SignatureGenes`������Signature�л��򼯣��ö���","�ָ�  
    * `corAnalysis_Method`������Է���������  
        * "pearson"  
        * "spearman"  
    * `calOrder`������Signatureʱ��׼����ѡ��  
        * standFirst��ÿ�������������������н��б�׼�����ټ���Signature Score  
        * calFirst���ȼ���Signature Score���ٽ��б�׼��  
        * noStand�������б�׼��  
    * `pickMethod`������Signature Score�ķ�����  
        * mean��ÿ��������Signature Genes��ƽ�����ֵ  
        * median��ÿ��������Signature Genes���м���ֵ  
    * `log`������ת��  
        * log2(n)  
        * log2(n+1)  
        * FALSE  
    * `exp`��ȡ��  
        * TRUE  
        * FALSE  


