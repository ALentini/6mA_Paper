## This function calculates the fraction of spliced reads in BAM files
# Reads aligned using STAR (2.6.0c) assign "N" in CIGAR if the read is spliced
# For paired-end reads, each mate is checked individually if any of the mates are spliced
# As splice junction locations may not be perfect, a small gap is allowed (in this case 5 bp)
# To get an expected distribution, all reads are resampled and checked for splicing

get.spliced <- function(x,spliceSites,paired=F,maxgap=5,n=1e5){
  require(GenomicAlignments)
  if(paired==T){
    message("Paired-end mode")
    find.splice <- function(aln){
      m1 <- grepl("N",cigar(aln@first),fixed = T)
      m2 <- grepl("N",cigar(aln@last),fixed = T)
      return(mean(m1+m2 > 0))
    }
    message("Filtering data")
    aln <- readGAlignmentPairs(x,strandMode = 0,param=ScanBamParam(flag=scanBamFlag(isProperPair = T)))
    spl <- aln[overlapsAny(aln@first,spliceSites,maxgap = maxgap) | overlapsAny(aln@last,spliceSites,maxgap = maxgap)]
    message("Calculating splice events")
    obs <- find.splice(spl)
    message("Bootstrapping")
    exp <- replicate(n,find.splice(sample(aln,size = length(spl))))
    list(obs=obs,exp=exp)
  }else{
    message("Single-end mode")
    find.splice <- function(aln){
      mean(grepl("N",cigar(aln),fixed = T))
    }
    message("Filtering data")
    aln <- readGAlignments(x)
    spl <- aln[overlapsAny(aln,spliceSites,maxgap = maxgap)]
    message("Calculating splice events")
    obs <- find.splice(spl)
    message("Bootstrapping")
    exp <- replicate(n,find.splice(sample(aln,size = length(spl))))
    list(obs=obs,exp=exp)
  }
}

## Exon-intron junctions
# This gets the exon-adjacent bases for introns of known transcript
# Non-exon-overlapping introns are retrived then disjoint against a 1 bp shrunk version of itself, resulting in 1 bp fragments flanking exons
# Visual representation:
# ========= --> intron
#  =======  --> intron-1
# =       = --> spl
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
intron <- keepStandardChromosomes(intronicParts(TxDb.Hsapiens.UCSC.hg38.knownGene),species = "Homo sapiens",pruning.mode = "coarse")
intron <- intron[width(intron)>1]
spl <- disjoin(c(intron,intron-1))
spl <- spl[width(spl)==1]

# Example: get.spliced(X_Aligned.sortedByCoord.out.bam,spliceSites=spl,paired=T)
