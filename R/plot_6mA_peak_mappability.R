options(stringsAsFactors = F)

####################
## Load peak data ##
####################

library(data.table)
meta <- fread("meta.tsv")[Sample %like% "DNA"]

library(rtracklayer)
pfls <- list.files("peaks_Bowtie2",pattern = ".narrowPeak.gz$",full.names = T)
pfls <- pfls[!grepl("Input",meta$Sample[match(gsub("_.*","",basename(pfls)),meta$Accession)])]

bed <- lapply(pfls,import.bed,extraCols=c(signalValue = "numeric", pValue = "numeric",  qValue = "numeric", peak = "integer")) 
bed <- lapply(bed,keepStandardChromosomes,species="Homo_sapiens",pruning.mode="coarse")
names(bed) <- gsub("_.*","",basename(pfls))

pfls_mc <- list.files("Peaks_Roadmap_Liftover",pattern=".narrowPeak.gz$",full.names = T)
bed_mc <- lapply(pfls_mc,import.bed,extraCols=c(signalValue = "numeric", pValue = "numeric",  qValue = "numeric", peak = "integer")) 
bed_mc <- lapply(bed_mc,keepStandardChromosomes,species="Homo_sapiens",pruning.mode="coarse")
names(bed_mc) <- gsub("_.*","",basename(pfls_mc))

##########################
## Mappability in peaks ##
##########################

# Calculate average mappability score over peaks
sum.score <- function(x,bw,fun=mean, ...){
  require(GenomicRanges)
  ol <- findOverlaps(x,bw)
  tapply(score(bw)[subjectHits(ol)],queryHits(ol),fun, ...)
}

# NOTE: these files are converted from doi:10.1093/nar/gky677, check the umap folder for more info!
# WARNING: The below code is resource-intensive!
umap50 <- import.bw("umap/k50.umap.bw")
umap100 <- import.bw("umap/k100.umap.bw")

# 6mA
map50 <- lapply(bed[names(bed) != "SRR6447606"],sum.score,bw=umap50,fun=mean, na.rm=T) # SRR6447606 = WGA
map50_avg <- sapply(map50,mean,na.rm=T)
map100 <- lapply(bed[names(bed) != "SRR6447606"],sum.score,bw=umap100,fun=mean, na.rm=T)
map100_avg <- sapply(map100,mean,na.rm=T)
# 5mC
map_mc50 <- lapply(bed_mc,sum.score,bw=umap50,fun=mean, na.rm=T)
map_mc50_avg <- sapply(map_mc50,mean,na.rm=T)
map_mc100 <- lapply(bed_mc,sum.score,bw=umap100,fun=mean, na.rm=T)
map_mc100_avg <- sapply(map_mc100,mean,na.rm=T)

# Generate shuffled background
library(BSgenome.Hsapiens.UCSC.hg38)
shuffle.genome <- function(x,BSgenome){
  require(GenomicRanges)
  si <- SeqinfoForBSGenome(BSgenome)
  lim <- seqlengths(si)[as.character(seqnames(x))] - width(x)
  pos <- round(runif(length(x),1,lim))
  GRanges(seqnames(x),IRanges(pos,width=width(x)))
}
set.seed(42) # reproducible regions
bg <- shuffle.genome(keepSeqlevels(Reduce(c,bed),paste0("chr",c(1:22,"X","Y")),pruning.mode = "coarse"),BSgenome.Hsapiens.UCSC.hg38)

map_bg50 <- replicate(100,get.score(sample(bg,round(mean(sapply(bed,length)))),bw = umap50,fun = mean,na.rm=T))
map_bg50_avg <- apply(map_bg50,2,mean,na.rm=T)
map_bg100 <- replicate(100,get.score(sample(bg,round(mean(sapply(bed,length)))),bw = umap100,fun = mean,na.rm=T))
map_bg100_avg <- apply(map_bg100,2,mean,na.rm=T)

##################
## Plot results ##
##################
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(cowplot)

# Average boxplot
df_map <- melt(list(k50=list("6mA"=map50_avg,"5mC"=map_mc50_avg,"Random"=map_bg50_avg),k100=list("6mA"=map100_avg,"5mC"=map_mc100_avg,"Random"=map_bg100_avg)))

ggplot(df_map,aes(y=value,x=L2)) + geom_boxplot(outlier.shape = NA) + geom_quasirandom(stroke=NA) + theme_bw() + ylab("Average mappability") + coord_cartesian(ylim=c(0,1)) + facet_wrap(~L1) + theme(axis.title.x = element_blank())

# Barplot
frac <- 0.5
df_map_pc <- melt(list(k50=list("mA"=lapply(map50,function(x) mean(x<frac,na.rm=T)),"mC"=lapply(map_mc50,function(x) mean(x<frac,na.rm=T))),k100=list("mA"=lapply(map100,function(x) mean(x<frac,na.rm=T)),"mC"=lapply(map_mc100,function(x) mean(x<0.5,na.rm=T)))))
df_map_bg <- data.frame(value=c(mean(apply(map_bg100,2,function(x) mean(x<frac,na.rm=T) )),mean(apply(map_bg50,2,function(x) mean(x<frac,na.rm=T) ))),L1=c("k100","k50"))
# with(df_map_pc,aggregate(value,list(L2,L1),mean)) # group average

ggplot(df_map_pc,aes(y=value,x=reorder(L3,value,sum),fill=L2)) + geom_col() + geom_hline(data=df_map_bg,aes(yintercept=value)) + facet_wrap(~L1) + coord_flip() + theme(legend.position = "top", axis.title.y = element_blank()) + ylab("Peak mappability < 50%")
