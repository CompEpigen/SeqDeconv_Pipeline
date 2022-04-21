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
dimension (k, number of cell types in this
    case).

    ##           [,1]       [,2]     [,3]      [,4]       [,5]      [,6]      [,7]
    ## [1,] 0.8994027 0.02053952 0.403561 0.5968592 0.90710214 0.5610518 0.3626625
    ## [2,] 0.1005973 0.97946048 0.596439 0.4031408 0.09289786 0.4389482 0.6373375
    ##           [,8]       [,9]     [,10]      [,11]        [,12]     [,13]     [,14]
    ## [1,] 0.4863976 0.97527458 0.1353662 0.95521844 0.0001114793 0.5704866 0.1373949
    ## [2,] 0.5136024 0.02472542 0.8646338 0.04478156 0.9998885207 0.4295134 0.8626051
    ##           [,15]     [,16]    [,17]     [,18]     [,19]      [,20]
    ## [1,] 0.09521017 0.8888249 0.833656 0.8541527 0.1197887 0.03918126
    ## [2,] 0.90478983 0.1111751 0.166344 0.1458473 0.8802113 0.96081874

    ## [1]  2 20

# **Generate pseudo-bulk samples**

Firstly, we calculated the ratio for the read sampling from pure
cell-type bulks according to the composition calculated above and a
total number of reads of the pseudo-bulk. Let’s assume we would like to
generate the first pseudo-bulk with 136,000,000 reads. The total number
of reads within each pure cell-type bulk is required (This can be
calculated using *samtools stats*).

    ## [1] 0.13705493 0.07341844

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
