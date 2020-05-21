# set directory to where the RNA and Clinical folders are
setwd("/Users/federicomatteo/Downloads/")
library(survival)

# read RNA file 
rna <- read.table('RNA/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
rna <- rna[-1,]

# and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
clinical <- t(read.table('Clinical/KIRC.merged_only_clinical_clin_format.txt', header=T, row.names=1, sep='\t'))

# first I remove genes whose expression is == 0 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))
  r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}
remove <- rem(rna)
rna <- rna[-remove,]

table(substr(colnames(rna),14,14))

n_index <- which(substr(colnames(rna),14,14) == '1')
t_index <- which(substr(colnames(rna),14,14) == '0')
