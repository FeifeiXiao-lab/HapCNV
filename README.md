# HapCNV
HapCNV: A Comprehensive Framework for CNV Detection in Haploid Single-Cell DNA Sequencing Data

# Author
Xuanxuan Yu, Fei Qin, Shiwei Liu, Noah Brown, Guoshuai Cai, Jennifer L. Guler, Feifei Xiao

# Description
Most CNV detection methods are specifically developed for diploid genomes without specific consideration of effects on haploid genomes. Single-cell or low-input sequencing data generally displays shallow and highly non-uniform read counts resulting from the whole genome amplification steps that introduce amplification biases. Additionally, reference samples or normal controls are used to provide baseline signals for defining copy number losses or gains. In traditional methods, references are usually pre-specified from cells that are assumed to be normal or disease-free. However, the use of pre-defined reference cells can bias results if common CNVs are present. Here, we present the development of a comprehensive statistical framework for data normalization and CNV detection in haploid single-cell or low-input DNA sequencing data called HapCNV. The predominant advancement is the construction of a novel genomic location specific pseudo-reference that selects unbiased references using a preliminary cell clustering method. Our approach effectively preserves common CNVs and is with easy applicability to diploids.

# Installation
HapCNV requires several R packages which requires manually installed: `modSaRa`,`FLCNA`,`GenomicRanges`, and `DNAcopy`. 

```
install.packages("devtools")
library(devtools)
install_github("FeifeiXiao-lab/HapCNV")
```

# Running HapCNV
## Examples

