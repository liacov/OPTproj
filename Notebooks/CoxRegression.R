# set directory to where the RNA and Clinical folders are
setwd("/Users/federicomatteo/Downloads/")
library(survival)

# read RNA file 
new <- read.table('RNA/KIRC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',nrows=20533, header=T,row.names=1,sep='\t')
# and take off first row cause we don't need it
new <- new[-1,]

# and read the Clinical file, in this case i transposed it to keep the clinical feature title as column name
clinical <- t(read.table('/Users/federicomatteo/Desktop/OPTproj/Data/Clinical/KIRC.merged_only_clinical_clin_format.txt', header=T, row.names=1, sep='\t'))
clinical = as.data.frame(clinical)
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
dim(rna)
write.table(rna, "/Users/federicomatteo/Desktop/OPTproj/Data/mydata.txt", sep=";")
n_index <- which(substr(colnames(new), 14, 14) == '1')
t_index <- which(substr(colnames(new), 14, 14) == '0')



sapply(rownames(rna), function(x) unlist(strsplit(x,'\\|'))[[1]])
colnames(rna) <- gsub('\\.','-',substr(colnames(rna),1,12))
View(rna3)
clinical$IDs <- toupper(clinical$patient.bcr_patient_barcode)
sum(df$IDs %in% colnames(rna3))

n_index <- which(substr(colnames(new),14,14) == '1')
t_index <- which(substr(colnames(new),14,14) == '0')
# save rna3 
rna[, n_index & t_index & df$IDs]
clinical$death_event <- ifelse(clinical$patient.vital_status == 'alive', 0,1)


ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}

clinical$death_days = death_collapsed

ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum(is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

clinical$followUp_days = fl_collapsed


clinical$new_death <- c()
for (i in 1:length(as.numeric(as.character(clinical$death_days)))){
    clinical$new_death[i] <- ifelse (is.na(as.numeric(as.character(clinical$death_days))[i]),
                                    as.numeric(as.character(clinical$followUp_days))[i],as.numeric(as.character(clinical$death_days))[i])
}

clinical$new_death

clinical$death_event

new_clinical = clinical[,c("new_death","death_event","IDs")]
new_clinical


plot(new_clinical$new_death, new_clinical$death_event)

write.table(new_clinical, "/Users/federicomatteo/Desktop/OPTproj/Data/SurvivalTimes.txt", sep=";")
