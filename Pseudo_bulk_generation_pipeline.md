Pseudo-bulk sample generation pipeline
================
Yunhee Jeong
20 April 2022

For the evaluation, we generated pseudo-bulk samples by randomly
sampling reads from pure cell-type samples based on given compositions.
We used Dirichlet distribution to simulate varied biological scenarios.

# **Cell-type composition simulation using Dirichlet distribution**

*MeDeCom* package provides a function called *generateExample*
simulating cell-type compositions based on Dirichlet distribution. This
function generate synthesised methylation array samples with dissected
components. Synthesised cell-type composition X bulk matrix is also
provided as A matrix. Please see the github page for the implementation
(<https://github.com/lutsik/MeDeCom/blob/master/R/utilities.R>).

*generateExample* requires three parameters: number of genomic features
(m, this doesn’t affect the A matrix so we just gave 10000), number of
profiles (n, number of bulks in this case) and number of hidden
dimension (k, number of cell types in this case).

``` r
library(MeDeCom, quietly = T)
example_res <- generateExample(m = 10000, n=20, k=2)
proportions <- example_res[["A"]] # Simulated proportions are saved as A matrix
head(proportions)
```

    ##           [,1]      [,2]      [,3]      [,4]         [,5]        [,6]      [,7]
    ## [1,] 0.8757685 0.8093913 0.3514253 0.7202147 9.999986e-01 0.008590047 0.9025135
    ## [2,] 0.1242315 0.1906087 0.6485747 0.2797853 1.387475e-06 0.991409953 0.0974865
    ##            [,8]       [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
    ## [1,] 0.96337204 0.01945097 0.7167934 0.3806838 0.8434721 0.2947496 0.1294949
    ## [2,] 0.03662796 0.98054903 0.2832066 0.6193162 0.1565279 0.7052504 0.8705051
    ##          [,15]      [,16]   [,17]    [,18]      [,19]     [,20]
    ## [1,] 0.6929803 0.95179319 0.71074 0.547832 0.03575168 0.5205605
    ## [2,] 0.3070197 0.04820681 0.28926 0.452168 0.96424832 0.4794395

``` r
dim(proportions)
```

    ## [1]  2 20

# **Generate pseudo-bulk samples**

Firstly, we calculated the ratio for the read sampling from pure
cell-type bulks according to the composition calculated above and a
total number of reads of the pseudo-bulk. Let’s assume we would like to
generate the first pseudo-bulk with 136,000,000 reads. The total number
of reads within each pure cell-type bulk is required (This can be
calculated using *samtools stats*).

``` r
proportion <- proportions[,1]
total_n_reads = 136000000
df_pure_cell_type <- data.frame(ctype = c("mL6-2", "mPv"),
                                n_reads = c(892479849, 186346020))
ctype_read_bulk <- (proportion * total_n_reads) / df_pure_cell_type$n_reads
ctype_read_bulk # Ratio for sampling reads from each pure cell-type bulk
```

    ## [1] 0.1334534 0.0906673

Read sampling can be done using *samtools*. You can give the ratio for
sampling with *-s*.

``` bash
samtools view -b -s 0.02936864 -o bulk_1_mL6-2.bam mL6-2.bam 
```

Once sampling of all cell types are done, the generated .bam files were
merged using *samtools merge*.

``` bash
samtools merge -r bulk_1.bam bulk_1_mL6-2.bam bulk_1_mPv.bam
samtools index bulk_1.bam
```
