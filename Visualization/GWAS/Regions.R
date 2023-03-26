library(plyr)
# note1 you can speed this up with the data.table package and changing read.table and write.table to fread and fwrite
# note2 if you already unzipped your info.gz files update this to .info
# note3 you can change the MAF and Rsq filter to whatever you think is appropriate 

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.20)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  da <- dat[,c("SNP","ALT_Frq","Rsq")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".info",sep = ""), row.names = F, quote = F, sep = "\t")
}

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.20)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
  da <- dat[,c("range")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}