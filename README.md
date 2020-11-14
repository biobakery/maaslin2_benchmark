Large-scale Benchmarking of Microbial Multivariable Association Methods
================

-   [Introduction](#introduction)
-   [Simulation Designs](#simulation-designs)
-   [Simulation Methods](#simulation-methods)
-   [Examples](#examples)
-   [Output](#output)
-   [Testing A New Method](#testing-a-new-method)
-   [Folder Descriptions ](#folder-descriptions)
-   [References](#references)
-   [Citation](#citation)

Introduction
------------
With the recent development of modern, massively parallel DNA sequencing technologies, there is a growing interest in developing robust and powerful statistical methods to uncover relationships between organisms and multivariate clinical covariates of their hosts (e.g. disease status, diet, and lifestyle variables, among others). Existing methods for microbe-metadata association have primarily focused on univariable (pairwise) association methods. While some of these methods can also perform multivariable association (e.g. [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [metagenomeSeq](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html), [limma-voom](https://bioconductor.org/packages/release/bioc/html/limma.html), and [MaAsLin](https://huttenhower.sph.harvard.edu/maaslin)), there is a surprising dearth of systematic studies aimed at multiple covariates and repeated measures, with no clear consensus on the appropriate method for scalable epidemiological studies. As microbiome data continue to accumulate, there is a pressing need for a rigorous comparison of these multivariable tools to help guide the researchers and practitioners alike. To this end, a large-scale benchmarking of microbial multivariable association tests is performed, essentially extending the differential abundance benchmarking of [McMurdie and Holmes (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531), [Weiss et al. (2017)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y), and [Hawinkel et al. (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28968702) to multivariable association analysis with or without repeated measures. 

The basic idea is to input a table of feature counts (e.g. taxa counts or gene/transcript/metabolite abundances on a per sample basis) and a table of metadata (such as sample-wise clinical covariates) and evaluate performance across different statistical models with varied normalization and transformation schemes. This is done in two stages. First, `Load_Simulator.R` is executed with desired parameters (details below) to generate synthetic datasets with user-defined spike-in information (both cross-sectional and repeated measures designs are supported). The null synthetic metagenomic counts are first generated using [sparseDOSSA](https://github.com/biobakery/sparseDOSSA), which are subsequently combined with the user-defined metadata information for spiking-in. Next, `Simulation_Main.R` script is run, which upon completion, creates an output file containing various performance metrics (e.g. False Discovery Rate (FDR), False Positive Rate (FPR), F1 Score, and Area Under the Curve (AUC), among others), based on replicated runs (typically hundreds of iterations). 

Simulation Designs
------------------

At the top level, we experiment with different metadata structures (**metadataType**). For each **metadataType**, we vary a combination of other parameters such as number of subjects (**nSubjects**), number of samples per subject (**nPerSubject**), number of microbes (**nMicrobes**), proportion of spiked-in microbes (**spikeMicrobes**), number of metadata (**nMetadata**), proportion of spiked-in metadata (**spikeMetadata**), and effect size (**effectSize**). The metadata structures are as follows:


| metadataType  | Type of Association               | Description                                 |
| ------------- |:---------------------------------:| -------------------------------------------:|
| UVA           | Univariable continuous            | One continuous metadata                     |
| UVB           | Univariable binary                | One dichotomous (binary) metadata           |
| MVA           | Multivariable mixed, independent  | 50% binary and 50% continuous, with no correlation between them |
| MVB           | Multivariable mixed, collinear    | 50% binary and 50% continuous, with AR(1) correlation structure (autocorrelation parameter = 0.5)|                           

For each of the above scenario (**metadataType**), we vary the metadata-agnostic parameters as follows:

* **nSubjects** `10 20 50 100 200`

* **nMicrobes** `100 200 500`

* **spikeMicrobes** `0.1`

* **effectSize** `1 2 5 10`

In addition, for multiple covariates, we vary the following metadata parameters (for univariable association these values are fixed at 1):

* **nMetadata** `5`

* **spikeMetadata** `0.2`

Similarly, for repeated mesaures, we vary the following parameter (fixed at 1 for cross-sectional designs):

* **nPerSubject** `5`


Simulation Methods
------------------

We use the following naming convention: `Model`.`Normalization`.`Transformation`, i.e. for each class of methods, we apply a set of models, that can take different input tables based on user-defined normalization and/or transformation methods. For example, `LM.TSS.LOG` means a linear model (LM) with TSS normalization and LOG transformation is desired. If no normalization/transformation is applied or a default pipeline is implemented, the corresponding entries are blank (e.g. metagenomeSeq, edgeR, DESeq2, etc.). Both fixed effects and random effects modeling are implemented for each method (when applicable). For GLM-based models, the logarithm of the normalized library size (with the exception of TSS for which original library size is retained) is included as an offset parameter. 

In summary, the following normalzation methods are tested:

```bash

CSS

RLE

TMM

TSS

CLR

```

In addition, the following transformation methods are tested:


```bash

AST

LOG

LOGIT

```

A complete list of tested methods are as follows:


```bash

ANCOM 

CPLM (Compound Poisson)

DESeq2

edgeR

limma

limma2 (limma with library size or scale factor as covariate)

limmaVOOM

LM

LM2 (LM with library size or scale factor as covariate)

metagenomeSeq

metagenomeSeq2

negbin (Negative Binomial)

Spearman

Wilcoxon

ZIB (Zero-inflated Beta)

ZINB (Zero-inflated Negative Binomial)

```

In addition to the above core scenarios, a subset of noncore scenarios are also evaluated to investigate the effect of (i) rarefaction, (ii) filtering, (iii) multiple testing adjustment methods, (iv) pseudocounts, (v) `nPerSubject`, (vi) `nMetadata`, (vii) zero-inflation, (viii) library size, (ix) multicollinearity, and (x) data-generating model.

Examples
--------

The sample codes are as follows:

```bash
Rscript Load_Simulator.R --metadataType UVB --nSubjects 100 --nPerSubject 5 --nMicrobes 200 --spikeMicrobes 0.1 --nMetadata 1 --spikeMetadata 1 --effectSize 20 
```

```bash
Rscript Simulation_Main.R --methodName LM.TSS --metadataType UVB --nSubjects 100 --nPerSubject 5 --nMicrobes 200 --spikeMicrobes 0.1 --nMetadata 1 --spikeMetadata 1 --effectSize 20 
```

Output
------

Running the above will produce an output file in the `Output` folder, containing various performance metrics for each test, which can be loaded in R and plotted and interpreted as desired. 

Testing A New Method
--------------------

To test a new method, one can add a function to run the test, and save it in the `Library` folder. As input, it should take two data.frame objects like above (with matched rows) containing features and metadata respectively and output a feature by metadata table of coefficients, p-values, and q-values. Examples can be found in the `Library` folder.


Folder descriptions 
-----------------

**Codes**
    - Two main simulation functions typically called from terminal with command-line arguments: `Load_Simulator.R` (for generating synthetic abundances, accompanied by spiked-in metadata) and `Simulation_Main.R` (univariable and multivariable tests).

**Library**
   - Utility functions and scripts for running individual methods.

**Input**
   - Folder to store all synthetic data sets generated by `Load_Simulator.R`.

**Output**
   - Folder to store output from `Simulation_Main.R`.

**Shell Scripts**
   - Examples on how to parallelize the scripts with bash to obtain all parameter combinations in a high performance computing cluster such as [Odyssey](https://www.rc.fas.harvard.edu/odyssey-3-the-next-generation/).

**Workflow Files**
- Examples on how to automate the above parallelization using [ANADAMA2](http://huttenhower.sph.harvard.edu/anadama2).

References
----------

McMurdie, P. J., and Holmes, S. (2014). [Waste not, want not: why rarefying microbiome data is inadmissible](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531). PLoS Computational Biology, 10(4):e1003531.

Weiss, S. et al. (2017). [Normalization and microbial differential abundance strategies depend upon data characteristics](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y). Microbiome, 5:27.

Hawinkel, S. et al. (2017). [A broken promise: microbiome differential abundance methods do not control the false discovery rate](https://www.ncbi.nlm.nih.gov/pubmed/28968702). Briefings in Bioinformatics, bbx104.


Citation
--------

Mallick, H. et al. (2020+). [Multivariable Association Discovery in Population-scale Meta-omics Studies](http://huttenhower.sph.harvard.edu/maaslin) (In Submission).

