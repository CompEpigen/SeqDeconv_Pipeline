MethylPurify cell-type deconvolution analysis
================
Yunhee Jeong
11 November 2021

MethylPurify detects differentially methylated regions (DMRs) and
estimated tumour purity from tumor methylome bulk samples. Reference
fasta file and the specification of CG island in .bed file are required
to run this software. Details of the usage is available here:
<https://pypi.org/project/MethylPurify/#description>

``` bash

ref="./bismark/hg19.fa"
cpg_bed="./ref/CGI_hg19_slop1000.bed"
out_dir="./bulk_1/"

bin_size=200
coverage=10
sampling=50

MethylPurify -o $out_dir -f bulk_1.bam -b $bin_size -c $coverage -s $sampling -i $cpg_bed -g $ref --species "hg19"
```

There are two output files. *MethylProfile.bed* contains detected DMRs
information and *alpha1.pred* shows the estimated tumour
    purity.

``` bash
cat MethylProfile.bed
```

    ## #chr start   end m1  m2  read_count  cytosine_count_in_CG    likelihold
    ## chr1 13146   13346   0.527089433533726   0.548177745041164   21  3   -24.131789023097646
    ## chr1 13349   13549   0.5189147967934834  0.9306057259456227  49  4   -66.78750195158416
    ## chr1 13556   13756   0.8809747309632764  0.8741802410686826  34  4   -51.37477183232292
    ## chr1 13759   13959   0.9967528407511198  0.9334398255438524  33  2   -7.942365830938414
    ## chr1 14544   14744   0.9928236929072323  0.9368540977469867  20  7   -23.946534988654783

``` bash
cat alpha1.pred 
```

    ## 0.255
    ## 6244
