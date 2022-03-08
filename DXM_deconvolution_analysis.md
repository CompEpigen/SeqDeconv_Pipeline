DXM cell-type deconvolution analysis
================
Yunhee Jeong
8 March 2022

ClubCpG is a computational tool extracting heterogeneous signals in WGBS
read-level.

# **Data preparation**

DXM requires a bed file containing chromosome, start, end, methylation
level columns, coverage and region\_id delimited by tab within a
designated set of regions. For our benchmarking, we used the proomter
and CGI regions extracted from mm10 and hg19 genomes each.

Firstly, we extracted methylation level at all CpG sites using
*MethylDackel* as follow. The user has to provide reference genome for
this step.

``` bash
bam="./bulk_1.bam"
f_out="./bulk_1"
MethylDackel extract ./ref/mm10.fa $bam -o $f_out
```

Then we loaded promoter regions from UCSC genome dataset and filtered
all CpG sites located on the loaded promoter regions.

``` r
library(GenomicRanges)
library(pbmcapply)
library(stringr)
library(dplyr)
library(data.table)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

COV <- 4 # coverage for filtering
NCORES <- 30 # Number of cores for multi-processing

# Promoters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promos <- promoters(txdb, upstream = 2000, downstream = 400)

# CpG island
# Downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/ and https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
cgi = fread(paste0(HOME_DIR, "Deconvolution_benchmark/DXM/cpgIslandExt_mm10.txt"), sep="\t", stringsAsFactors = F, header = F, showProgress = T)
colnames(cgi) <- c("no", "chr","start", "end", "ncpg", "V6", "V7", "V8", "V9", "V10", "V11")
cgi = makeGRangesFromDataFrame(cgi)

bedgraph <- "./bulk_1_CpG.bedGraph"
df_bulk <- fread(bedgraph, sep="\t", skip = 1, stringsAsFactors = F, header = F, nThread = NCORES, data.table = F)
colnames(df_bulk) <- c("chr", "start", "end", "methyl_level", "m", "n")

# Calculate coverage
df_bulk$cov = df_bulk$m + df_bulk$n
gr_bulk = makeGRangesFromDataFrame(df_bulk, keep.extra.columns = T)

# Calculate overlap
df_final = mclapply(1:length(promos), function(i){
    #print(promos[i,])
    df_intersect = GenomicRanges::findOverlaps(promos, gr_bulk)
    if(length(df_intersect)){
        df_intersect = gr_bulk[df_intersect@to,]
        df_intersect$region = paste0("region_", as.character(i))
        return(as.data.frame(df_intersect))
    }
    
}, mc.cores = NCORES)

# remove regions without any detected CpGs
df_final = df_final[lengths(df_final) != 0]

# Filter regions with the given coverage 
df_final = df_final[df_final$cov > COV, ]


# Save the result
bed_file = "bulk_1_DXM_promoter.bed"
df_final$methyl_level= (df_final$methyl_level)/100 # divide by 100 to make the number in the range between 0 and 1
write.table(df_final[, c("seqnames", "start", "end", "region", "methyl_level", "cov")], 
            file=bed_file, sep="\t", col.names = T, quote = F, row.names = F)
```

You can load the CGI regions after downloading the CpG Island file from
UCSC database.

``` r
# CpG island
# Downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/ and https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
cgi = fread("./cpgIslandExt.txt"), sep="\t", stringsAsFactors = F, header = F, showProgress = T)
colnames(cgi) <- c("no", "chr","start", "end", "ncpg", "V6", "V7", "V8", "V9", "V10", "V11") # Give column names only to the information we will use
cgi = makeGRangesFromDataFrame(cgi)
```

# **Running DXM**

DXM can be installed following the instruction on the github
(<a href="url">https://github.com/CompEpigen/dxm</a> ). It provides a
command line function to estimate cell-type composition from the
generated *.bed*. The explanation of each parameter can be found on the
github as well.

``` bash
dxm_estimateFracs -i bulk_1_promoter_CGI.bed -k 3 -o testPrevalence 
```

Then, you can get the estimation result in the file
*testPrevalence\_solvedPrevalences.txt*.

``` bash
0.10530123327336183
0.26727870626758243
0.6274200604590557
```
