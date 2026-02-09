# USERM  <a href="https://github.com/xiangmingcai/USERM"><img src="./man/figures/USERMlogo.png" align="right" width="100"/></a>

Unmixing Spread Estimation based on Residual Model

<!-- badges: start -->
[![R-CMD-check](https://github.com/xiangmingcai/USERM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xiangmingcai/USERM/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/xiangmingcai/USERM/branch/main/graph/badge.svg)](https://app.codecov.io/gh/xiangmingcai/USERM)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://xiangmingcai.github.io/USERM)

<!-- badges: end -->


## 🔍 Introduction

The USERM provides an out-of-box tool to apply the residual model approach [^1], which characterizes and predicts the spread of unmixed spectral flow cytometry data, which arises from instrumental noise or deviations between actual cellular emission and the average fluorescence signatures.

The USERM also supports computing various matrixes tools for panel design, including the Coef Matrix, the Hotspot Matrix, and the Similarity Matrix, and many others.

## 🔧 Updates
**2026-01-28**: 

1. Release version 1.0.1
2. Update instructions

## 💻 Installation

You can install the development version of USERM from [GitHub](https://github.com/xiangmingcai/USERM) with:

``` r
#install dependency with BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("flowCore", "ComplexHeatmap", "circlize"))

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("xiangmingcai/USERM")
```

## 🖊 Instruction

The following instructions and all reference can be found here: 📖 [Documentation website](https://xiangmingcai.github.io/USERM)

1. Basic use of the USERM

2. Using custom single-color control FCS file

3. Using USERM to interpret observed spread

## ⚙️ Contribute Single-color control FCS files to the USERM

The capabililty of the USERM greatly rely on the built-in FCS. However, there are many fluorescence not included in the current version. Also, there are many other spectral flow cytometry instruments that are missing in the current USERM. The current built-in SCC are all human PBMC or beads. That means we are also short in SCC acquired from other tissues or species.

It would be of great help if you would contribute your SCC to the USERM package. We will acknowledge your contribution in the "Resource" column of the querySig() return list. If you are interested in contributing your SCC, please contact us directly via e-mail: 

Dr. Juan J. Garcia Vallejo: jj.garciavallejo@amsterdamumc.nl

Xiangming Cai: x.cai@amsterdamumc.nl or xiangming_cai@126.com

## 🧮 Feedback
Any feedback is welcomed! You may post any issue on the [Issue](https://github.com/xiangmingcai/USERM/issues) of the USERM repository.


## 📚 Citation

If you use this package in your research, please cite our paper:
```
Xiangming Cai, Sara Garcia-Garcia, Leo Kuhnen, Michaela Gianniou, Juan J. Garcia Vallejo. Unmixing Spread Estimation Based on Residual Model in Spectral Flow Cytometry. 
bioRxiv 2026.01.27.701929; doi: https://doi.org/10.64898/2026.01.27.701929
```
