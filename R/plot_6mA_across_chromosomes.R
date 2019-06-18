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

#####################################
## Peak density across chromosomes ##
#####################################

library(BSgenome.Hsapiens.UCSC.hg38)
si <- keepStandardChromosomes(SeqinfoForBSGenome(BSgenome.Hsapiens.UCSC.hg38),"Homo_sapiens")

bins <- unlist(tileGenome(si,tilewidth = 1e6,cut.last.tile.in.chrom = T)) # 1Mbp bins

cnt <- sapply(bed,countOverlaps,query=bins)
cnt_rel <- apply(cnt,2,function(x) unlist(tapply(x,seqnames(bins), function(y) y/sum(y) )) )

cnt_mc <- sapply(bed_mc,countOverlaps,query=bins)
cnt_mc_rel <- apply(cnt_mc,2,function(x) unlist(tapply(x,seqnames(bins), function(y) y/sum(y) )) )

# NOTE: RMSK bed file from UCSC: http://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=232602379_xRfRuviMPL722frTYr1gpAEtsJ2N&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&hgta_outputType=primaryTable&hgta_outFileName=hg38_rmsk.bed
rmsk <- fread("zcat hg38_rmsk.bed.gz")
l1_gr <- keepStandardChromosomes(with(rmsk[repFamily=="L1"],GRanges(genoName,IRanges(genoStart,genoEnd),strand=strand)),"Homo_sapiens","coarse")
cnt_l1 <- unlist(tapply(countOverlaps(bins,l1_gr),seqnames(bins),function(x) x/sum(x) ))

library(ggplot2)
library(reshape2)
chr <- factor(as.character(seqnames(bins)),levels=paste0("chr",c(1:22,"X","Y","M")))
df_melt <- melt(cnt_rel)
df_melt$sample <- meta$Sample[match(df_melt$Var2,meta$Accession)]
df_melt <- rbind(df_melt,cbind(melt(cnt_mc_rel),sample="DNA_mC"))
df_melt$x <- start(bins)
df_melt$chr <- chr
df_melt <- df_melt[-grep("chrM",df_melt$Var1),]

# NOTE: Assembly gaps from UCSC: http://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=232602379_xRfRuviMPL722frTYr1gpAEtsJ2N&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=gap&hgta_table=0&hgta_regionType=genome&hgta_outputType=primaryTable&hgta_outFileName=gap.txt
gaps <- fread("zcat gap.txt.gz")
df_gap <- gaps[chrom %in% standardChromosomes(BSgenome.Hsapiens.UCSC.hg38),c(2:4,8)]
colnames(df_gap) <- c("chr","x1","x2","class")
df_gap$chr <- factor(df_gap$chr,levels=paste0("chr",c(1:22,"X","Y")))

df_l1 <- melt(cnt_l1)
df_l1$x <- start(bins)
df_l1$chr <- chr
df_l1 <- df_l1[-grep("chrM",df_l1$chr),]

# Plot genome-wide
ggplot(df_melt[-grep("WGA",df_melt$sample),],aes(y=value,x=x/1e6,col=sample,fill=sample)) + geom_smooth() + facet_wrap(~chr, scales = "free") + labs(y="Fraction peaks", x="Mbp") + theme(legend.position = "top") + geom_line(data=df_l1,inherit.aes = F,col="black",aes(y=value,x=x/1e6)) + geom_rect(data=df_gap,alpha=0.5, inherit.aes = F,aes(ymin=0,ymax=0.01,xmin=x1/1e6,xmax=x2/1e6,fill=class))

# Plot ideogram
library(Gviz)
ideo <- IdeogramTrack(genome = "hg38")
lapply(seq_along(si),function(i) plotTracks(ideo,chromosome=seqnames(si)[i],from=1,to=seqlengths(si)[i]) )


#############################
## Genome-wide correlation ##
#############################

# Remove bin 3103 (chrM)
mean(cor(cnt_rel[-3103,],cnt_l1[-3103]))
mean(cor(cnt_mc_rel[-3103,],cnt_l1[-3103]))

library(ggplot2)
library(ggbeeswarm)
library(cowplot) # for theme
# Plot correlation
ggplot(data.frame(y=c(cor(cnt_rel[-3103,],cnt_l1[-3103]),cor(cnt_mc_rel[-3103,],cnt_l1[-3103])),x=rep(c("6mA","5mC"),c(ncol(cnt_rel),ncol(cnt_mc_rel)))),aes(y=y,x=x,col=x)) + geom_quasirandom(stroke=NA) + geom_boxplot(outlier.shape = NA,col="black",fill=NA) + coord_cartesian(ylim=c(0,1)) + ylab("Genome-wide LINE1 correlation (r)") + theme(axis.title.x=element_blank()) + guides(col=F)
