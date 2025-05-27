# AiES: Axon Integlity Evaluation System

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD* is a free and open-source R package for (high-throughput) axonal integrity analysis of micrographs, such as murine dorsal root ganglia explant cultures.
AiES offers tools to segment neurites, extract quantitative features, generate SVM models, and quantify axon integrity.
It is released under the BSD 3-Clause License.

## Installation

```r
# Install from CRAN (when available)
# install.packages("AiES")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("BreezyCave/AiES")
```

## Overview

AiES provides an integrated workflow for axon image analysis, including feature extraction, machine learning classification, and quantitative assessment.
The main functions are:

1. `axDistmap()`: Extracts morphological features from axon images (TIFF image files).
2. `axSvm()`: Builds an SVM classifier utilizing extracted features.
3. `axQnt()`: Quantifies Axon Integrity Index (AII) and Degeneration Index (DI) for new data using a trained SVM.

For detailed usage instructions, please refer to the vignettes included in the package.

## Typical Workflow

### 1. Feature Extraction from Images

```r
library(AiES)

# Specify sample image folders included in the package
img_dir1 <- system.file("extdata", "Degenerate_Images", package = "AiES")
img_dir2 <- system.file("extdata", "Intact_Images", package = "AiES")

# Extract features from multiple folders (recursive search)
axDistmap(
  folder_paths = c(img_dir1, img_dir2),
  subBack = 30,
  resizeW = 900,
  output_path = tempdir()
)

```
### Exporting All Features

```r
axDistmap(
  folder_paths = img_dir1,
  allFeatures = TRUE,
  output_path = tempdir()
)
```

This function processes all TIFF files in the selected directories, generating for each image:
- A distance map image
- A binary image (optional)
- A text file containing feature data

### 2. Building an SVM Classifier

```r
Degenerate_dir <- system.file("extdata", "Degenerate_txt", package = "AiES")
Intact_dir <- system.file("extdata", "Intact_txt", package = "AiES")

axSvm(
  degenerate_path = Degenerate_dir,
  intact_path = Intact_dir,
  output_data_path = tempdir(),
  output_model_path = tempdir(),
  nCst = 3, nGmm = 0.1, nCrss = 5
)

```
This function generates:
- An extracted data file for machine learning
- A trained SVM model file

### 3. Quantitative Analysis
```r
img_dir1 <- system.file("extdata", "Degenerate_Images", package = "AiES")
img_dir2 <- system.file("extdata", "Intact_Images", package = "AiES")

result <- axQnt(
  input_dir = c(img_dir1, img_dir2),
  svm_model_path = system.file("extdata", "svm_example_model.svm", package = "AiES"),
  output_dir = tempdir()
)

# Visualize the results
library(ggplot2)
ggplot(result, aes(x = AxonIntegrityIndex)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  ggtitle("AII Distribution")

```
This function generates:
- A summary CSV file containing the Axon Integrity Index and Degeneration Index for each image
- (Optional) Single image prediction data files

### License
AiES is licensed under the BSD 3-Clause License.
See the LICENSE file for details.

### Citation

If you use AiES in your research, please cite:

Tokunaga S, Funakoshim M, Araki T, et al. AiES: Axon integrity Evaluation System. National Center of Neurology and Psychiatry, 2025.


## Contact

For bug reports or feature requests, please use the Issues page on GitHub.
