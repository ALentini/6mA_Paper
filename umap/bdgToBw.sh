## Converts Umap bedgraph files to bigwig
# fetchChromSizes and bedGraphToBigWig are from UCSC: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
# Umap tracks are from doi:10.1093/nar/gky677 (https://bismap.hoffmanlab.org/)
# Umap k50 multi-mappability track: https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz
# Umap k100 multi-mappability track: https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bedgraph.gz

fetchChromSizes hg38 > hg38.genome #| awk '$1 ~ /chr.{1,2}$/' | sort -k1,1 -V -s

for file in *.umap.bedgraph.gz
do
	zcat $file | tail -n +2 | LC_COLLATE=C sort -k1,1 -k2,2n > ${file%.*}.sorted.bedgraph
	bedGraphToBigWig ${file%.*}.sorted.bedgraph hg38.genome ${file%%.*}.umap.bw
	rm ${file%.*}.sorted.bedgraph
done
