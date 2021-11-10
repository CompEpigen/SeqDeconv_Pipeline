csmFinder\_coMethy\_deconvolution\_analysis
================
Yunhee Jeong
10 November 2021

## Data Preparation

To use csmFinder and coMethy, appropriate format of data has to be
generated via *bismark\_methylation\_extractor*. It generates a cytosine
report from given bam file and reference fasta
file.

``` bash
bismark_methylation_extractor --comprehensive -o CpG_context_bulk_1.txt.gz --gzip --cytosine_report --genome_folder ./hg19/bismark/ bulk_1.bam
```

## csmFinder

csmFinder extracts putative cell-subset specific DNA methylation (pCSM)
loci from methylomes saved in the bismark cytosine report file. Details
of each step are in
<a href="url"><https://github.com/Gavin-Yinld/csmFinder>\<\>.

``` r
library(bedtoolsr)
library(csmFinder)
if(!requireNamespace("stringr")) BiocManager::install("stringr")
library(stringr)

f_cpg = "hg19_plus_cpg_cytosine.txt"
f_bulk <- "CpG_context_bulk_1.txt.gz"

segment <- bismark2segment(files = f_bulk, CpG_file = f_cpg, tmp_folder="./tmp")
candidate <- find_candidate(segment, depth = 10)
pcsm_segment <- csmFinder(candidate,distance=0.3,pval=0.05)
pcsm_loci <- merge_segment(pcsm_segment)
```

``` r
head(segment)
```

    ##                                             V1                    V2
    ## 1     chr1:10000005_10000072_10000273_10000320               1111:1;
    ## 2     chr1:10000072_10000273_10000320_10000327               1111:1;
    ## 3 chr1:100000827_100000912_100000960_100001014        0101:1;1101:1;
    ## 4 chr1:100000912_100000960_100001014_100001022 0111:1;1001:1;1011:3;
    ## 5 chr1:100000960_100001014_100001022_100001073 0000:1;0010:3;1110:1;
    ## 6 chr1:100001014_100001022_100001073_100001199        0000:1;1111:2;

``` r
head(candidate)
```

    ##                                              V1
    ## 44 chr1:100008901_100008946_100008951_100008991
    ## 67 chr1:100015941_100015992_100016002_100016012
    ## 69 chr1:100015992_100016002_100016012_100016047
    ## 70 chr1:100016002_100016012_100016047_100016054
    ## 71 chr1:100016012_100016047_100016054_100016107
    ## 72 chr1:100016047_100016054_100016107_100016110
    ##                                      V2
    ## 44               0000:10;0111:1;1111:2;
    ## 67        0000:18;1100:1;1110:1;1111:1;
    ## 69               0000:16;1101:1;1111:2;
    ## 70 0000:16;0001:1;1000:1;1011:1;1111:3;
    ## 71               0000:12;0010:1;1111:3;
    ## 72               0000:13;0100:1;1111:7;

``` r
head(pcsm_segment)
```

    ##                                               V1
    ## 71  chr1:100016012_100016047_100016054_100016107
    ## 72  chr1:100016047_100016054_100016107_100016110
    ## 74  chr1:100016107_100016110_100016121_100016157
    ## 122 chr1:100025156_100025200_100025239_100025289
    ## 300 chr1:100066190_100066206_100066318_100066341
    ## 318 chr1:100070552_100070587_100070598_100070618
    ##                                       V2         d         pval
    ## 71                0000:12;0010:1;1111:3; 0.9807692 1.036890e-08
    ## 72                0000:13;0100:1;1111:7; 0.9821429 4.769877e-10
    ## 74  0000:11;0001:1;1101:1;1110:1;1111:6; 0.9166667 3.935477e-02
    ## 122         0000:8;0001:1;0111:1;1111:5; 0.9305556 2.845974e-02
    ## 300         0000:8;0001:1;1000:1;1111:4; 0.9500000 2.482661e-05
    ## 318               0000:1;1000:13;1111:3; 0.7678571 4.769877e-10

``` r
head(pcsm_loci)
```

    ##     V1     V2     V3
    ## 1 chr1 540839 540877
    ## 2 chr1 564501 564691
    ## 3 chr1 566287 566305
    ## 4 chr1 566340 566466
    ## 5 chr1 566618 566714
    ## 6 chr1 567206 567358

## coMethy

For coMethy, we collected all pCSM loci from all bulk samples and
converted the sequencing methylation data into an array format of
samples by available pCSM
    loci.

``` r
print(head(bulk_mat))
```

    ##                       bulk_1    bulk_2    bulk_3    bulk_4    bulk_5    bulk_6
    ## chr1_20129_20234   0.9583333 0.8757576 0.9555556 0.9649123 0.8675214 0.8273810
    ## chr1_62096_62156   0.7134810 0.8270771 0.7542088 0.7597222 0.8455270 0.8465909
    ## chr1_63627_63683   0.7500000 0.9137255 0.7500000 0.7500000 0.8863636 0.8873016
    ## chr1_76896_76995   1.0000000 0.9226190 1.0000000 0.9166667 0.9000000 0.9166667
    ## chr1_102827_102995 0.5306638 0.8709273 0.4853098 0.5501089 0.8203396 0.8486111
    ## chr1_128786_128846 0.3333333 0.7592593 0.5000000 0.3000000 0.6888889 0.6878307
    ##                       bulk_7    bulk_8    bulk_9   bulk_10   bulk_11   bulk_12
    ## chr1_20129_20234   0.8273810 0.9052288 0.8681917 0.8995098 0.9103314 0.8379085
    ## chr1_62096_62156   0.8465909 0.8411462 0.8335116 0.8690476 0.8447060 0.8315508
    ## chr1_63627_63683   0.8952991 0.8611111 0.8717949 0.9722222 0.9761905 0.9050926
    ## chr1_76896_76995   0.9166667 0.9000000 0.9166667 0.9500000 0.9583333 0.9166667
    ## chr1_102827_102995 0.8639077 0.7536659 0.8150183 0.7943355 0.7659938 0.8632278
    ## chr1_128786_128846 0.6878307 0.6888889 0.6878307 0.9393939 0.7797619 0.7012987
    ##                      bulk_13   bulk_14   bulk_15   bulk_16   bulk_17   bulk_18
    ## chr1_20129_20234   0.8668047 1.0000000 0.8634921 0.9203431 0.9285714 0.9152778
    ## chr1_62096_62156   0.8514212 0.8095517 0.8860119 0.8478706 0.7261905 0.8271505
    ## chr1_63627_63683   0.9218501 0.7797619 0.9351852 0.9761905 0.8888889 0.9743590
    ## chr1_76896_76995   0.9131944 0.9375000 0.9000000 0.9392857 1.0000000 0.9305556
    ## chr1_102827_102995 0.8634080 0.5914352 0.9095960 0.7727920 0.5765432 0.7618899
    ## chr1_128786_128846 0.7595960 0.6111111 0.7555556 0.8769841 0.4539683 0.8560606
    ##                      bulk_19   bulk_20
    ## chr1_20129_20234   0.9047619 0.8534572
    ## chr1_62096_62156   0.7960834 0.8882258
    ## chr1_63627_63683   0.9722222 1.0000000
    ## chr1_76896_76995   0.9187500 0.9583333
    ## chr1_102827_102995 0.7661836 0.8732875
    ## chr1_128786_128846 0.8424242 0.9444444

``` r
print(dim(bulk_mat))
```

    ## [1] 454605     20

coMethy package is comprised of two steps. Firstly, it finds hypo-, mid-
and hyper-methylatio groups through k-means clustering. Secondly, from
these clusters, it makes co-methylation module and extract eigen-loci
from respective module. Again, details are explained in
<a href="url"><https://github.com/Gavin-Yinld/coMethy>\<\>.

``` r
kmeans_clust <- co_methylation_step1(bulk_mat)
module <- co_methylation_step2(data=bulk_mat, 
                                 kmeans_result = kmeans_clust, 
                                 softPower_list = sp_list, plot=T) ## 26,28,30 for 2 cell types

eigen_loci <- extract_eigen(methy_data = module$profile, all_label = module$module_id, number_of_eig=n_eigens, plot=T)
```

Then, the final cell-type prporortion estimation is done by MeDeCom
using the
eigen-loci.

``` r
medecom.result<-runMeDeCom(as.matrix(eigen_loci$methy_prof), Ks = 2:3, 10^(-5:-2), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=CORES)
proportion <- MeDeCom::getProportions(medecom.result, K=5)
```

``` r
proportion
```

    ##      bulk_1 bulk_2 bulk_3 bulk_4    bulk_5    bulk_6   bulk_7    bulk_8
    ## LMC1      0      1      0      0 0.3629549 0.6880197 0.651378 0.1760309
    ## LMC2      1      0      1      1 0.6370451 0.3119803 0.348622 0.8239691
    ##        bulk_9   bulk_10   bulk_11 bulk_12 bulk_13    bulk_14 bulk_15   bulk_16
    ## LMC1 0.377046 0.5127846 0.3842423       1       1 0.02966781       1 0.6591912
    ## LMC2 0.622954 0.4872154 0.6157577       0       0 0.97033219       0 0.3408088
    ##        bulk_17  bulk_18   bulk_19    bulk_20
    ## LMC1 0.0382916 0.575407 0.4981166 0.96267137
    ## LMC2 0.9617084 0.424593 0.5018834 0.03732863
