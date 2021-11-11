MeDeCom cell-type deconvolution analysis
================

MeDeCom is an array-based cell-type deconvolution method based on
non-negative matrix factorisation (NMF). In order to apply this method
to our sequencing data, we converted the data into the array shape and
chose 2000 CpGs with the highest methylation level varience.

# ****Example data****

``` r
dim(sample_mat)
```

    ## [1] 20000     3

``` r
head(sample_mat)
```

    ##      bulk_1_CpG bulk_2_CpG bulk_3_CpG
    ## [1,]          1          0          1
    ## [2,]          1          0          1
    ## [3,]          0          1          0
    ## [4,]          0          1          0
    ## [5,]          0          1          0
    ## [6,]          0          1          0

# ****MeDeCom analysis****

You can run MeDeCom using *runMeDeCom* from MeDeCom package. We gave
smaller numbers for some parameters in this document due to running
time, but you can find original parameter set up in our
paper.

``` r
medecom_res <- runMeDeCom(as.matrix(sample_mat), Ks=2:3, lambdas=10^(-5:-2), NCORES=8, NFOLDS = 5, ITERMAX = 10, NINIT = 30)
```

    ## [Main:] checking inputs
    ## [Main:] preparing data
    ## [Main:] preparing jobs
    ## [Main:] 168 factorization runs in total
    ## [Main:] finished all jobs. Creating the object

*getProportions* function culculate cell-type composition from the
result of *runMeDeCom*. MeDeCom creates a result with latent DNA
methylation components (LMCs). You can consider this LMC as a cell type
in this analysis.

``` r
proportions <- getProportions(medecom_res)
colnames(proportions) <- colnames(sample_mat)
proportions
```

    ##      bulk_1_CpG   bulk_2_CpG bulk_3_CpG
    ## LMC1          1 9.666168e-14          1
    ## LMC2          0 1.000000e+00          0
