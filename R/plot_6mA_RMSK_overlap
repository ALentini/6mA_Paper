options(stringsAsFactors=F)

#####################
## Load peak files ##
#####################

meta <- read.delim("meta.tsv")

library(rtracklayer)
pfls <- list.files("peaks_Bowtie2",pattern = ".narrowPeak.gz$",full.names = T)
pfls <- pfls[!grepl("Input",meta$Sample[match(gsub("_.*","",basename(pfls)),meta$Accession)])]

bed <- lapply(pfls,import.bed,extraCols=c(signalValue = "numeric", pValue = "numeric",  qValue = "numeric", peak = "integer")) 
bed <- lapply(bed,keepStandardChromosomes,species="Homo_sapiens",pruning.mode="coarse")
names(bed) <- gsub("_.*","",basename(pfls))

# NOTE: RMSK bed file from UCSC: http://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=232602379_xRfRuviMPL722frTYr1gpAEtsJ2N&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&hgta_outputType=primaryTable&hgta_outFileName=hg38_rmsk.bed
rmsk <- keepStandardChromosomes(with(fread("zcat hg38_rmsk.bed.gz"),GRanges(genoName,IRanges(genoStart,genoEnd),strand=strand)),"Homo_sapiens","coarse")

# hg38 as a genomic reference
library(BSgenome.Hsapiens.UCSC.hg38)
si <- keepStandardChromosomes(SeqinfoForBSGenome(BSgenome.Hsapiens.UCSC.hg38),"Homo_sapiens")
rc_len <- tapply(rmsk,rmsk$repClass,function(x){sum(as.numeric(width(x)))})

#############################
## Plot total RMSK overlap ##
#############################

library(reshape2)
library(ggplot2)
library(ggbeeswarm)
df_sum <- data.frame(Var=gsub("_.*","",basename(pfls)),value=sapply(bed,function(x) mean(x %over% rmsk) ))
df_sum$Condition <- apply(na.omit(meta[match(df_sum$Var,meta$Accession),c("Sample","Study")]),1,paste0,collapse="_")
df_sum <- rbind(df_sum,data.frame(Var="Genomic",value=sum(as.numeric(width(rmsk)))/sum(as.numeric(seqlengths(si))),Condition="Genomic"))

# Plot
ggplot(df_sum,aes(x=Condition,y=value,fill=Condition)) + stat_summary(fun.y = "mean",geom="bar",position=position_dodge(0.9)) + geom_beeswarm(dodge.width = 0.9,stroke=NA) + ylab("Fraction peaks in RMSK") + coord_cartesian(ylim=c(0,1)) + theme(axis.title.x = element_blank()) + theme_bw()

###################################
## Plot overlap per repeat class ##
###################################

rc_cnt <- do.call(rbind,tapply(rmsk,rmsk$repClass,function(y){sapply(bed,function(x) mean(x %over% y) )}))
colnames(rc_cnt) <- gsub("_peaks.*","",basename(pfls))
rc_cnt <- cbind(rc_cnt,Genomic=rc_len/sum(as.numeric(seqlengths(si))))
rc_red <- rc_cnt[table(rmsk$repClass) > 1e4,] # Keep only major repeat classes
df <- melt(t(rc_red))
df$Condition <- "Genomic"; 
df$Condition[!grepl("Genomic",df$Var1)] <- apply(na.omit(meta[match(df$Var1,meta$Accession),c("Sample","Study")]),1,paste0,collapse="_")
df <- df[!grepl("Input",df$Condition),]

# Plot
ggplot(df,aes(x=Var2,y=value,fill=Condition)) + stat_summary(fun.y = "mean",geom="bar",position=position_dodge(0.9)) + geom_beeswarm(dodge.width = 0.9,stroke=NA) + geom_text(data=df[df$value > 0.3 & df$Var2 == "LINE",],aes(label=Var1,x=Var2,y=value)) + ylab("Fraction in repeat") + theme(axis.title.x = element_blank()) + theme_bw()

##########################################
## Plot paired WGA and Native DIP peaks ##
##########################################

# Plot
ggplot(df[df$Var1 %in% c("SRR6447606","SRR6447607"),],aes(x=Var2,y=value,fill=Var1)) + geom_col(position="dodge") + theme(legend.position = "top") + theme_bw()
# Correlation
cor.test(rc_red[,"SRR6447606"],rc_red[,"SRR6447607"])

# Fraction overlapping peaks
tapply(rmsk,rmsk$repClass,function(z){
  x <- bed[["SRR6447607"]]
  y <- bed[["SRR6447606"]]
  xidx <- x %over% z
  yidx <- y %over% z
  return(c(mean(x[xidx] %over% y[yidx]),mean(y[yidx] %over% x[xidx])))
})
