# set directory to where the RNA and Clinical folders are
setwd("/Users/federicomatteo/Downloads/")
library(survival)

# read RNA file 
new <- read.table('RNA/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
new <- new[-1,]

# and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
clinical <- t(read.table('Clinical/KIRC.merged_only_clinical_clin_format.txt', header=T, row.names=1, sep='\t'))
df = as.data.frame(clinical)
names(df)
# first I remove genes whose expression is <= 1 in more than 50% of the samples:
rem <- function(x){
  x <- as.matrix(x)
  x <- t(apply(x, 1, as.numeric))
  r <- as.numeric(apply(x, 1, function(i) sum(i <= 275.735))) 
  remove <- which(r > dim(x)[2]*0.5)
  return(remove)
}

rna = new
remove <- rem(rna)
rna <- rna[-remove,]
dim(rna)
table(substr(colnames(rna),14,14))

n_index <- which(substr(colnames(new), 14, 14) == '1')
t_index <- which(substr(colnames(new), 14, 14) == '0')


sapply(rownames(rna), function(x) unlist(strsplit(x,'\\|'))[[1]])
colnames(rna) <- gsub('\\.','-',substr(colnames(rna),1,12))
View(rna3)
df$IDs <- toupper(df$patient.bcr_patient_barcode)
sum(df$IDs %in% colnames(rna3))

n_index <- which(substr(colnames(new),14,14) == '1')
t_index <- which(substr(colnames(new),14,14) == '0')
# save rna3 
rna[, n_index & t_index & df$IDs]
all_clin$death_event <- ifelse(df$patient.vital_status == 'alive', 0,1)

