# ����Signature��GO����ע�ͺ�ͨ·֮�������Է���  

## ���  

�˴����ǽ�Ϊ�û��Ѻ��͵Ľű�������**ĳһ��������**��**���ܡ�ͨ·**֮���**����Է���**�������н϶�ɹ�ѡ��Ĳ������ҿ����Զ��������Ͻ�������Է�����  
PS��Ҳ��������**��������**�빦�ܡ�ͨ·֮�������Է�����  
���⻹�ṩ�˽�����ͳһת��ΪSymbol�Լ��������ж��ֲ�����ϵ�������ű���  

## �����밲װ  

ֱ�����ز��ڱ���ʹ��  

## ����  

+ GeneSignature-Pathways_Correlation.Rmd������Է���  
+ GSEA_GeneSets_ID2Symbol.Rmd��GeneID����ת��ΪSymbol  
+ GSEA_GeneSets_Normalize2Symbol.Rmd��GeneName����ת��ΪSymbol�Ľű������Ƽ�
+ loop.R�����ֲ�����ϵ���������Է���
+ .\lib\GeneSignature-Pathways_Correlation_Function.R��GeneSignature-Pathways_Correlation.Rmd�������Զ��庯��  

## ������ʹ��  

### һ��GeneSignature-Pathways_Correlation.Rmd  

#### ��һ������

�ô���ͨ��ʹ�û���������Ϣ������ĳһ��Gene Signature��GO����ע�͡�Pathwayͨ·֮�������Է�����Gene Signature��һ���������ѧ����Ļ��򼯺ϣ�ͨ��������򼯵ı������Signature�Ĵ�֣��������Signature�������еı�������  
��ps��֧�����뵥�����򣬼���û���������ܡ�ͨ·֮��������

#### ������ʹ��

1. ����  
**ע��**�����нű���Ҫ������ԭʼ���ݷ����ڹ���Ŀ¼�µ�data�ļ�����  
    1. ��������ݣ�����TCGA��CGGA�����ݣ������������������  
    **����Ҫ��**����Ϊ������Ϊ��������Ϊcsv��ʽ  
        + TCGA���ݣ����Դ�TCGA�����������ݲ�ͨ��script��TCGA_Rawdata_Merge-Expression.Rmd������������ʽ���  
        + CGGA���ݣ����Դ�CGGA�����������ݲ�ͨ��Project��[CGGA_RawData_Merge-Combat](https://github.com/SuiYinG2000/CGGA_RawData_Merge-Combat)����������ʽ���  
        + �������������  
    2. Gene���ƶ�ӦSymbolע���ļ�  
    �Ƽ���NCBI GENE������[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    3. GO����ע�ͺ�ͨ·����  
    **����Ҫ��**����Ϊgmt��ʽ  
    �Ƽ���GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)��վ������**Gene ID**��ɵĻ����ļ�  
    ��ԭ��������Gene Symbols���ļ������򼯻�������վ��չʾ��������ͬʱ���ֻ��򲢲���Symbol��  
    **ע��**��
        + ��Ҫͬʱ����**GO����ע��**��**Pathway**�Ļ��򼯣������ڲ��ֶ��޸Ĵ���Ӷ�������������һ��������ԡ�  
        + ����ͨ��**GSEA_GeneSets_ID2Symbol.Rmd**����**GSEA_GeneSets_Normalize2Symbol.Rmd**������ͳһת��ΪSymbol  
    4. �����Զ��庯����  
.\lib\GeneSignature-Pathways_Correlation_Function.R��GeneSignature-Pathways_Correlation.Rmd  

2. ����  
**����**����·����Ϣ��Ϊ����·���������ýű������ڹ���Ŀ¼��˫���򿪣�����ѡ��ʹ�����·��
    + `workdir`������Ŀ¼  
    + `libdir`�������Զ��庯����  
    + `expression_file`�����������  
      + "\\data\\XXX.csv"  
    + `GOGeneSet`��GO����
      + "\\data\\XXX.gmt"
    + `PWGeneSet`��Pathway����  
      + "\\data\\XXX.gmt"
    + `GeneSet_Name_Normalized`�����򼯻��������Ƿ�ͨ��GSEA_GeneSets_ID2Symbol.Rmd��GSEA_GeneSets_Normalize2Symbol.Rmd����ת��  
      + TRUE  
      + FALSE
    + `GeneNameAnno`��Gene���ƶ�ӦSymbolע���ļ�  
      + "\\data\\XXX.csv"
    + `SignatureGenes_Name`������Signature�����ƻ򵥸���������  
      + "m6A_Signature_Genes"  
      + "STK26"  
    + `SignatureGenes`������Signature�л��򼯣��ö���","�ָ�  
    + `corAnalysis_Method`������Է���������  
      + "pearson"  
      + "spearman"  
    + `calOrder`������Signatureʱ��׼����ѡ��  
      + "standFirst"��ÿ�������������������н��б�׼�����ټ���Signature Score  
      + "calFirst"���ȼ���Signature Score���ٽ��б�׼��  
      + "noStand"�������б�׼��  
    + `pickMethod`������Signature Score�ķ�����  
      + "mean"��ÿ��������Signature Genes��ƽ�����ֵ  
      + "median"��ÿ��������Signature Genes���м���ֵ  
    + `log`������ת��  
      + "log2(n)"  
      + "log2(n+1)"  
      + FALSE  
    + `exp`��ȡ��  
      + TRUE  
      + FALSE  

### ����GSEA_GeneSets_ID2Symbol.Rmd

#### ��һ������

�ô������ڽ���GESA MSigDB��վ�����ص�**Gene IDs**��ɵĻ��򼯣�ͳһת����**Gene Symbol**����ʽ  

#### ������ʹ��

1. ����  
**ע��**�����нű���Ҫ������ԭʼ���ݷ����ڹ���Ŀ¼�µ�data�ļ�����  
    1. Gene���ƶ�ӦSymbolע���ļ�  
    �Ƽ���NCBI GENE������[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    2. GO����ע�ͺ�ͨ·��GeneID����  
    **����Ҫ��**����Ϊgmt��ʽ  
    ��GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)��վ������**Gene ID**��ɵĻ����ļ�  
    **ע��**��
        + ��Ҫͬʱ����**GO����ע��**��**Pathway**�Ļ��򼯣������ڲ��ֶ��޸Ĵ���Ӷ�������������һ��������ԡ�  
2. ����  
    + `workdir`������Ŀ¼  
    + `libdir`�������Զ��庯����  
    + `resultsdir`������ļ��У����������
    + `picsdir`��ͼƬ�ļ��У����������
    + `GOGeneset`��GO����
      + "\\data\\XXX.gmt"
    + `PWGeneSet`��Pathway����  
      + "\\data\\XXX.gmt"
    + `GeneNameAnno`��Gene���ƶ�ӦSymbolע���ļ�  
      + "\\data\\XXX.csv"

### ����GSEA_GeneSets_Normalize2Symbol.Rmd

#### ��һ������

�ô������ڽ���GESA MSigDB��վ�����ص�**Gene Symbols**��ɵĻ��򼯣�ͳһת����**Gene Symbol**����ʽ  

#### ������ʹ��

1. ����  
**ע��**�����нű���Ҫ������ԭʼ���ݷ����ڹ���Ŀ¼�µ�data�ļ�����  
    1. Gene���ƶ�ӦSymbolע���ļ�  
    �Ƽ���NCBI GENE������[Homo_sapiens.gene_info.csv](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz)  
    2. GO����ע�ͺ�ͨ·��GeneID����  
    **����Ҫ��**����Ϊgmt��ʽ  
    ��GSEA [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb)��վ������**Gene Symbols**��ɵĻ����ļ�  
    **ע��**��
        + ��Ҫͬʱ����**GO����ע��**��**Pathway**�Ļ��򼯣������ڲ��ֶ��޸Ĵ���Ӷ�������������һ��������ԡ�  
2. ����  
    + `workdir`������Ŀ¼  
    + `libdir`�������Զ��庯����  
    + `resultsdir`������ļ��У����������
    + `picsdir`��ͼƬ�ļ��У����������
    + `GOGeneset`��GO����
      + "\\data\\XXX.gmt"
    + `PWGeneSet`��Pathway����  
      + "\\data\\XXX.gmt"
    + `GeneNameAnno`��Gene���ƶ�ӦSymbolע���ļ�  
      + "\\data\\XXX.csv"

### �ġ�loop.R

#### ��һ������

�ýű�Ϊ����Է�����������ű����������벻ͬ������Ͻ����������㡣

#### ������ʹ��

1. ����
    + `rmd`��GeneSignature-Pathways_Correlation.Rmd����λ�ú͸�rmd�ļ�
    + `workdir`������Ŀ¼
    + `libdir`�������Զ��庯����Ŀ¼
    + `expr_file`�����������
      + "\\data\\XXX.csv"
    + `GOGeneset`��GO����
      + "\\data\\XXX.gmt"
    + `PWGeneSet`��Pathway����  
      + "\\data\\XXX.gmt"
    + `GeneSet_Name_Normalized`�����򼯻��������Ƿ�ͨ��GSEA_GeneSets_ID2Symbol.Rmd��GSEA_GeneSets_Normalize2Symbol.Rmd����ת��  
      + TRUE  
      + FALSE
    + `GeneNameAnno`��Gene���ƶ�ӦSymbolע���ļ�  
      + "\\data\\XXX.csv"
    + `SignatureGenes_Name`������Signature�����ƻ򵥸���������  
      + "m6A_Signature_Genes"  
      + "STK26"  
    + `SignatureGenes`������Signature�л��򼯣��ö���","�ָ�  

    + `corAnalysis_Methods`������Է��������������ַ���/������
      + "pearson"  
      + "spearman"  
    + `calOrders`������Signatureʱ��׼����ѡ�񣨵����ַ���/������
      + "standFirst"��ÿ�������������������н��б�׼�����ټ���Signature Score  
      + "calFirst"���ȼ���Signature Score���ٽ��б�׼��  
      + "noStand"�������б�׼��  
    + `pickMethods`������Signature Score�ķ����������ַ���/������
      + "mean"��ÿ��������Signature Genes��ƽ�����ֵ  
      + "median"��ÿ��������Signature Genes���м���ֵ  
    + `exps`��ȡ�ݣ������ַ���/������
      + TRUE  
      + FALSE  
    + `logs`������ת���������ַ���/������
      + "log2(n)"  
      + "log2(n+1)"  
      + FALSE  

## ����

���Project������˶ʿ�����е�һ�ݹ��������ڲο�������û�����Ƶ����̽��ܣ�����һ�е�һ�ж��Ǵ��㿪ʼ������������TCGA�ͣ�CGGA���������ؼ�Ԥ����ȥ����ЧӦ����������ô����Signature�Ĵ�֣�median/exp���������һƪ��ʦ�Ƹ���������ʹ�õķ�������log(n+1)��GEPIA2��ʹ�õķ�����Ԥ�Ƚ��б�׼��������ȡƽ��ֵ��ФФ��ʦ�����ĺܰ���ָ�����Դ������Ų��ɣ��������������Է�����ʹ�á���ʱ�����¾��������Ĺ��̣�hh����Ҳû���ر�ļ�����������������������ôһ��Protocal��  

������������ƪREADME�ĵ�Ϊֹ�����׵�ͼҲû���ܹ���ȫ�ĸ��ֳ����ġ�ϲ�˵��ǣ�ʹ������Ϊ��Ϊ��ȷ��һ�ײ�����Ϻ󣬾�TCGA��CGGA���ݼ������ǿ������ԵĹ�����ͨ·����һ�£���ǿ������Ե�overlap�����ǿ���⣬CGGA�����ռ���õ��ĸ�����Ե�����Ҳ���١���������Ҫ������������������ϣ�CGGA��������̫��ÿ�μ���һ��������Ͼ�Ҫ3��4h����orz����  

��֮���ܹ��������������һ��Project���Ķ��ٻ�����һ���Ժ��еģ���˵����������ֵ�ʱ����ǵ��Լ��Ĵ��������Ժܸߣ����Ǻ��Ϳ�������һ�����迴�Ŵ����ں�̨�Լ������ܣ��ͺ������Լ������˵�С����һ����  

���������Ȥ������Ƚ�ˮ�Ĵ��벢���һЩ������������һ�ǳ���л���:)��  

��������Peace out.
