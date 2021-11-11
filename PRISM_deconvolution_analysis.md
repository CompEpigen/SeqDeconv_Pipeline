PRISM cell-type deconvolution analysis
================

PRISM infers epigenetically heterogeneous subclones from bisulfite
sequencing data saved as BAM file. According to authors, this tool has
been verified only with RRBS data. Detailed usage is described in
<a href="url">https://github.com/dohlee/prism</a>

# ****Run PRISM****

Prism firstly extracts epiloci where a group of reads are mapped in a
short genomic region. Default settings for read coverage 20 (d=20) and
minimum number of CpGs 4 (c = 4) are given in *prism extract command*.
*bulk\_1.bam* and *bulk\_1.met* file names are given input and ouput
files each.

``` bash
prism extract -i bulk_1.bam -o bulk_1.met -u -v
```

We corrected errorneous methylation pattern in produced *bulk\_1.met*
file through *prism preprocess* command which is based on hidden Markov
model. The result is saved in *bulk\_1\_corrected.met* file.

``` bash
prism preprocess -i bulk_1.met -o bulk_1_corrected.met
```

prism deconvolute command infers the subclonal composition of the
sample. Simply give methylation pattern-corrected epiloci file.

Finally, we conducted deconvolution on *bulk\_1\_corrected.met* file
using *prism deconvolute* command. Maximum number of cluster is set as
10 (-m option) and *bulk\_1.result* is the final result file.

``` bash
prism deconvolute -i bulk_1_corrected.met -o bulk_1.result -m 10 
```

# ****Cacluate cell-type proprtion from the PRSIM result****

Since *PRISM* only infers subclone of each epilocus, we calculated
proportions of each subclone as the fianl cell-type composition.

This is *bulk\_1.result* file readed in R as a data
    frame:

``` r
head(dt_res)
```

    ##                                                      epilocus cluster subclone
    ## 1 chr1:240426862;chr1:240426875;chr1:240426912;chr1:240426944       1        1
    ## 2     chr2:79940406;chr2:79940428;chr2:79940464;chr2:79940486       1        1
    ## 3     chr7:71041799;chr7:71041804;chr7:71041834;chr7:71041841       1        1
    ## 4 chr1:242570514;chr1:242570540;chr1:242570545;chr1:242570574       1        1
    ## 5 chr15:53716758;chr15:53716776;chr15:53716788;chr15:53716811       1        1
    ## 6     chr8:70602381;chr8:70602390;chr8:70602399;chr8:70602419       1        1
    ##   depths fingerprint_counts fingerprint_fractions\n    
    ## 1     19                  9                   0.4736842
    ## 2     21                  2                   0.0952381
    ## 3     20                  7                   0.3500000
    ## 4     18                  8                   0.4444444
    ## 5     21                  7                   0.3333333
    ## 6     23                 11                   0.4782609

We calculated ratio of each subclone and assigned the ratio as a
proportion of each cell type as below:

``` r
n_subclones = length(unique(dt_res$subclone))
composition <- unlist(lapply(unique(dt_res$subclone), function(i){nrow(dt_res[dt_res$subclone == i,])/nrow(dt_res)}))

if(sum(composition) == 0){
  composition <- unlist(lapply(-1:0, function(i){nrow(dt_res[dt_res$subclone == i,])/nrow(dt_res)}))
}

names(composition) <- lapply(1:n_subclones, function(x){paste0("ctype_", as.character(x))})
print(composition)
```

    ##   ctype_1   ctype_2 
    ## 0.8438215 0.1561785
