# AiES: Axon Integrity Evaluation System

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

**AiES** is a free and open-source R package for (high-throughput) axonal integrity analysis of micrographs, such as murine dorsal root ganglia explant cultures.
AiES offers tools to segment neurites, extract quantitative features, generate SVM models, and quantify axon integrity.
It is released under the BSD 3-Clause License.

## System requirements for Linux

AiES depends on EBImage, which needs the FFTW, JPEG, PNG and TIFF libraries. Before installing EBImage or AiES on Debian/Ubuntu, please run:

```
sudo apt-get update
sudo apt-get install -y libfftw3-dev libjpeg-dev libpng-dev libtiff-dev
```

Other distributions provide the same libraries under different package names. See the [EBImage page](https://bioconductor.org/packages/EBImage) for details.

## Installation

AiES requires R (>= 3.6.0). EBImage is distributed by Bioconductor, so it must be installed first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")

# Install from CRAN
install.packages("AiES")

# Or install the development version from GitHub
# install.packages("devtools")
# devtools::install_github("BreezyCave/AiES")
```

## Overview

AiES provides an integrated workflow for axon image analysis, including feature extraction, machine learning classification, and quantitative assessment.
The main functions are:

1. `axDistmap()`: Extracts morphological features from axon images (TIFF image files).
2. `axSvm()`: Builds an SVM classifier using the extracted features.
3. `axQnt()`: Quantifies the Axon Integrity Index (AII) and Degeneration Index (DI) using one of the two built-in models or a model trained with `axSvm()`.

AiES is designed for phase-contrast images in which axons appear darker than the background.
Images are resampled to the scale at which the built-in models were trained (0.791746 um per pixel), so the pixel size of the input should be known.
See `?AiES` for supported and not recommended image formats, and the vignette (`vignette("Introduction_to_AiES")`) for detailed usage.

## Quick start

```r
library(AiES)

# Sample images included in the package (900 px wide, 0.791746 um per pixel)
img_dir1 <- system.file("extdata", "Degenerate_Images", package = "AiES")
img_dir2 <- system.file("extdata", "Intact_Images", package = "AiES")

# Quantify the Axon Integrity Index (AII) and Degeneration Index (DI)
result <- axQnt(
  model = "g0.1_APA",                  # or "g0.2_ACC"
  input_dirs = c(img_dir1, img_dir2),  # multiple folders supported
  pixelSize = 0.791746,                # pixel size of the input images (um)
  output_dir = tempdir()
)
head(result)
```

To extract features (`axDistmap()`) and to train a classifier on your own data (`axSvm()`), see the vignette.

## Built-in models

| Model | Features | gamma | cost |
|---|---|---|---|
| `g0.1_APA` (recommended) | Ras, ENT, Per, SVA, SV2 | 0.0871 | 1.189 |
| `g0.2_ACC` | CON, ecc, mjx, cir, Ras, Ran | 0.207 | 1.110 |

## Upgrading from version 0.99.x

Feature names, the arguments of `axDistmap()`, and the output of `axQnt()` have changed, and AII and DI of the same images differ from those of version 0.99.x.
See the vignette (section "Upgrading from AiES 0.99.x") and `NEWS.md` for details.

### License

AiES is licensed under the BSD 3-Clause License.
See the LICENSE file for details.

### Notice regarding dependencies

This package makes use of the following R packages, each with its own license:

- EBImage: LGPL (>= 2.1)
- data.table: MPL-2.0
- e1071: GPL-2 | GPL-3
- fs: MIT
- tiff: GPL-2 | GPL-3
- rstudioapi: MIT

Users must comply with the license terms of these dependencies when redistributing or modifying this package.
See each package's CRAN/Bioconductor page for full license texts.

### Citation

If you use AiES in your research, please cite:

Tokunaga S, Funakoshi M, Araki T (2026). AiES: Axon Integrity Evaluation System for Microscopy Images. R package version 1.0.0. https://CRAN.R-project.org/package=AiES. https://doi.org/10.32614/CRAN.package.AiES

The same information is available in R with `citation("AiES")`.

## Contact

For bug reports or feature requests, please use the Issues page on GitHub.
