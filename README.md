# No evidence for DNA N6-methyladenine in mammals
### Douvlataniotis K*, Bensberg M*, Lentini A*, Gylemo B & Nestor CE
Availability: [Science Advances Paper](https://advances.sciencemag.org/content/6/12/eaay3335)

In this paper we show that all commonly used methods to detect DNA N6-methyladenine have inherent technical biases resulting in false positive detection.

This repository contains scripts and processed data used for the analyses in the paper.

### Pre-processing
GRCh38 genome idices were build using `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`.

For standard analysis, data was aligned using `bowtie2 [-N 1 --local]` (2.3.4.3).

Peaks were called against input using `macs2 [-g hs]` (2.1.2).

For analysis of read splicing, `STAR [--outSAMstrandField intronMotif --outSAMattribues NH HI AS nM XS]` (2.6.0c) was used.

### Data:
- `./meta.tsv` contains xIP-seq sample annotations for data included in the paper.
- `./Peaks_Bowtie2/` contains gzipped peaks called from Bowtie2 aligned data.
- `./Peaks_Roadmap_Liftover/` contains gzipped Roadmap 5mC DIP-seq peaks in GRCh38 coordinates.

### umap:
- `umap/bdgToBw.sh` was used to convert genome mappability tracks (doi:10.1093/nar/gky677) to bigWig format.

### R:
- `R/count_spliced_reads.R` was used to calculate fraction of spliced reads over splice sites.
- `R/plot_6mA_RMSK_overlap.R` was used to plot fraction peaks over repetitive features.
- `R/plot_6mA_across_chromosomes.R` was used to plot peak density across chromosomes.
- `R/plot_6mA_peak_mappability.R` was used to plot peak mappability (see `umap/bdgToBw.sh`).
