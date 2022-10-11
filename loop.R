# rm(list = ls())
# .rs.restartR()

rmd <- "E:\\GithubRepository\\Signature-Correlation\\GeneSignature-Pathways_Correlation.Rmd"

workdir <- "E:\\PostGraduate\\Botao\\Code_Verification\\SignatureGenes_GO-PW_Correlation"
libdir <-"E:\\GithubRepository\\Signature-Correlation\\lib"

# expr_files <- c("TCGA_Brain-GBM_AllSample_Gene-Name_Expression_fpkm_unstranded_2022-09-06.csv",
#                 "TCGA_Brain-GBM_AllSample_Gene-Name_Expression_tpm_unstranded_2022-09-27")
expr_file <- "\\data\\TCGA_Brain-GBM_AllSample_Gene-Name_Expression_fpkm_unstranded_2022-09-06.csv"

GOGeneset <- "\\data\\goGeneSet_normalized.gmt"
PWGeneset <- "\\data\\pathwayGeneSet_normalized.gmt"
GeneSet_Name_Normalized <- TRUE
GeneNameAnno <- "\\data\\Homo_sapiens.gene_info.csv"

SignatureGenes_Name <- "m6A_Signature_Genes"
SignatureGenes <- "METTL3,METTL14,WTAP,RBM15,ZC3H13,KIAA1429,METTL16,FTO,ALKBH5"

corAnalysis_Methods <- c("pearson", "spearman")
calOrders <- c("calFirst", "standFirst")
pickMethods <- c("median", "mean")
logs <- c(TRUE, FALSE)
exps <- c(TRUE, FALSE)


for (corMethod in corAnalysis_Methods){
  for (calOrder in calOrders){
    for (pickMethod in pickMethods){
      for (log in logs){
        for (exp in exps){
          if (!(log==T & exp==T)){
            tryCatch({
              params <- list(workdir = workdir,
                             libdir = libdir,
                             expression_file = expr_file,
                             GOGeneset = GOGeneset,
                             PWGeneset = PWGeneset,
                             GeneSet_Name_Normalized = GeneSet_Name_Normalized,
                             GeneNameAnno = GeneNameAnno,
                             SignatureGenes_Name = SignatureGenes_Name,
                             SignatureGenes = SignatureGenes,
                             corAnalysis_Method = corMethod,
                             calOrder = calOrder,
                             pickMethod = pickMethod,
                             log = log,
                             exp = exp)
              rmarkdown::render(rmd, params=params, envir=new.env())
              gc()
            }, warning = function(w){
              gc()
              next
            }, error = function(e){
              gc()
              next
            })
          }
        }
      }
    }
  }
}