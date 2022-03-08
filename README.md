# Sequencing-based DNA methylome deconvolution pipelines

![DeconvProcedure!](/figure/deconvolution_procedure.png "DeconvProcedure")


## Introduction 

Analysis of cell-type heterogeneity based on cell mixture (bulk) omic profiles is an active area of research, known as cell-type deconvolution. Sequencing-based DNA methylation data such as whole genome bisulfite sequencing (WGBS), in particular, has high capacity for cell-type deconvolution through leveraging read-level information. We thoroughly evaluated five previously published decovolution methods for sequencing data: Bayesian epiallele detection (BED), PRISM, csmFinder + coMethy, ClubCpg and MethylPurify, together with two array-based methods, MeDeCom and Houseman, as a comparison group. This github page describes the pipeline of each method including all the details (e.g. input data format, how to analyse the output results).

## Our evaluation results

**NOW OUR WORK IS ON BIORXIV!** Please find the results and our manuscript here [https://doi.org/10.1101/2021.11.29.470374
](https://doi.org/10.1101/2021.11.29.470374)



## References

- __BED:__ [James E Barrett, Andrew Feber, Javier Herrero, Miljana Tanic, Gareth A Wilson, Charles Swanton, and Stephan
Beck. Quantification of tumour evolution and heterogeneity via bayesian epiallele detection. BMC bioinformatics,
18(1):1–10, 2017.](https://doi.org/10.1186/s12859-017-1753-2)

- __ClubCpG:__ [C Anthony Scott, Jack D Duryea, Harry MacKay, Maria S Baker, Eleonora Laritsky, Chathura J Gunasekara,
Cristian Coarfa, and Robert A Waterland. Identification of cell type-specific methylation signals in bulk whole
genome bisulfite sequencing data. Genome biology, 21(1):1–23, 2020.](https://doi.org/10.1186/s13059-020-02065-5)

- __csmFinder+coMethy:__ [Liduo Yin, Yanting Luo, Xiguang Xu, Shiyu Wen, Xiaowei Wu, Xuemei Lu, and Hehuang Xie. Virtual methylome
dissection facilitated by single-cell analyses. Epigenetics & chromatin, 12(1):1–13, 2019.](https://doi.org/10.1186/s13072-019-0310-9)

- __DXM:__ [Fong, Jerry, et al. "Determining subpopulation methylation profiles from bisulfite sequencing data of heterogeneous samples using DXM." Nucleic acids research 49.16 (2021): e93-e93.](https://doi.org/10.1093/nar/gkab516)

- __MethylPurify:__ [Xiaoqi Zheng, Qian Zhao, Hua-Jun Wu, Wei Li, Haiyun Wang, Clifford A Meyer, Qian Alvin Qin, Han Xu,
Chongzhi Zang, Peng Jiang, et al. Methylpurify: tumor purity deconvolution and differential methylation detection
from single tumor dna methylomes. Genome biology, 15(7):1–13, 2014.](https://doi.org/10.1186/s13059-014-0419-x)

- __PRISM:__ [Dohoon Lee, Sangseon Lee, and Sun Kim. Prism: methylation pattern-based, reference-free inference of subclonal
makeup. Bioinformatics, 35(14):i520–i529, 2019.](https://doi.org/10.1186/s13059-014-0419-x)

- __Houseman's method:__ [Eugene Andres Houseman, William P Accomando, Devin C Koestler, Brock C Christensen, Carmen J Marsit,
Heather H Nelson, John K Wiencke, and Karl T Kelsey. Dna methylation arrays as surrogate measures of cell
mixture distribution. BMC bioinformatics, 13(1):86, 2012.](https://doi.org/10.1186/1471-2105-13-86)

- __MeDeCom:__ [Pavlo Lutsik, Martin Slawski, Gilles Gasparoni, Nikita Vedeneev, Matthias Hein, and Jörn Walter. Medecom:
discovery and quantification of latent components of heterogeneous methylomes. Genome biology, 18(1):1–20,
2017.](https://doi.org/10.1186/s13059-017-1182-6)
