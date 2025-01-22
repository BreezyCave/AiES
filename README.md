# Axon Integrity Index Calculation Package

This package provides a set of functions for calculating the Axon Integrity Index and Degeneration Index from axon images.

## Installation

```r
devtools::install_github("YourGitHubUsername/YourPackageName")
```

## Usage

This package includes three main functions:

1. `axDistmap()`: Extracts features from TIFF image files.
2. `axSvm()`: Creates an SVM model using the extracted features.
3. `axQnt()`: Calculates the Axon Integrity Index and Degeneration Index using the created SVM model.

For detailed usage instructions, please refer to the vignettes included in the package.

## Requirements

- R version 4.0.0 or higher
- The following R packages:
  - stringr
  - dplyr
  - EBImage
  - e1071
  - ggpubr

## License

This project is licensed under the BSD 3-Clause License. See the LICENSE file for details.

## Contact

For bug reports or feature requests, please use the Issues page on GitHub.
