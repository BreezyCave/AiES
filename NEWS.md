# AiES 1.0.1

This is a major update. Feature names, arguments and outputs have changed, so please read "Breaking changes" before upgrading from version 0.99.x.

## Breaking changes

* `axDistmap()`: the argument `allFeatures` is replaced by `featureSet = "standard"` or `"all"`.
* All features are now exported under short names (e.g. `s.perimeter` is `Per`, `h.ent.s1` is `ENT`, `s.area` is `Area`). See "Feature names" in `?axDistmap`. Text files written by version 0.99.6 (original EBImage names) are converted automatically by `axSvm()` and `axQnt()`.
* `axDistmap()`: the standard output now has 11 columns (Area, Per, Ras, Ran, ecc, mjx, ENT, SVA, SV2, CON, cir) instead of 7. Feature tables exported with the earlier default settings lack features required by the built-in models; please re-run `axDistmap()` on the images.
* The example model `svm_example_model.svm` is discontinued and replaced by two built-in models, `g0.1_APA` (recommended) and `g0.2_ACC`. AII and DI therefore differ from those of version 0.99.x.
* Images are now resampled to a physical scale (0.791746 um per pixel, the scale of the training data) instead of a fixed width of 900 px. The fixed width is used only when the pixel size is unknown, and a warning is issued.
* `axSvm()`: the default features and hyperparameters have changed (previously six fixed features, C = 3, gamma = 0.1). With `model = "custom"` and without `nCst` / `nGmm`, C = 1 and gamma = 0.1 are used.
* `axSvm()` now samples the training data by default: up to 12 images per group and up to 1000 objects per image, randomly selected with `seed = 1`. Training results therefore differ from version 0.99.6. To use all data, set `maxNumImg = NULL` and `maxNumObj = NULL`.
* `axSvm()`: `degenerate_path` and `intact_path` accept one folder or a character vector of .txt files. Several folders, or folders mixed with files, now stop with an error.
* `axQnt()`: the output is written to one subfolder per input folder (`<folder>_<date>_output_files`), and `FileName` in the summary is the path relative to the input folder.

## New features

* `axDistmap()` and `axQnt()`: new arguments `pixelSize`, `targetPixelSize`, `autoParam` and `colorMode`.
  * The pixel size is read from the `pixelSize` argument (recommended) or from the file metadata (OME-XML, TIFF resolution tags, ImageJ calibration).
  * Pixel-unit parameters are rescaled to the working scale when `autoParam = TRUE`.
  * `colorMode = "rgb"` averages the three channels of color-camera images.
  * Data with fewer significant bits than the TIFF container (e.g. 12-bit data in a 16-bit file) are rescaled automatically when the bit depth is recorded in OME-XML. A warning is issued for very dark images.
  * Binary (black-and-white) input images are detected and skipped.
* `axDistmap()` writes `<folder>_<date>_Parameters.csv` (one row per image). When analyzing images, `axQnt()` writes the same kind of file for each input folder.
* `axSvm()`: new arguments `model` (`"g0.1_APA"`, `"g0.2_ACC"` or `"custom"`), `features`, `targetPixelSize`, `maxNumImg`, `maxNumObj` and `seed`. The random number state of the R session is restored on exit. The scale of the training data is recorded with the model when `*_Parameters.csv` is found next to the input data.
* `axSvm()`: columns that are not present in every input file are reported and not used.
* `axQnt()`: new arguments `model` and `distMapImg`. A model produced by `axSvm()` is selected with `svm_model_path`, and the required predictors are read from the model itself. Models trained with version 0.99.x can still be used.
* `axQnt()`: when images were analyzed at a scale different from that of the model, one warning at the end of the run reports the number of images and the scale ratio.
* New help page `?AiES` describing supported, not recommended and unsupported images.

## Bug fixes

* `axQnt()`: fixed a missing binarization inversion in the image processing for images (`imprtImg = TRUE`), which made the segmentation differ from that of `axDistmap()`.
* `axQnt()`: the default of `svm_model_path` is now `NULL` (it was `tempdir()`).
* `axQnt()`: objects stored in the model file (e.g. a saved workspace) could silently overwrite the arguments of the function. The file is now loaded into a separate environment.
* `axQnt()`: a folder, a missing file or several paths given as `svm_model_path` now stop with a clear message. The model selection dialog accepts one file only.
* `axQnt()`: no summary file is written when no file could be processed.
* `axDistmap()`: `imgType = "jpg"` now works (converted to the `"jpeg"` accepted by `EBImage::writeImage()`).
* `axSvm()`: fixed errors when the folder selection dialog was canceled (R >= 4.3), when a vector of file paths was given (R >= 4.2), and in the save dialog of tcltk.
* `axSvm()`: input files with different columns (e.g. exported with different `featureSet` settings) no longer cause an `rbind()` error or silently dropped rows.

## Other changes

* Sample images: a second image, with the pixel size recorded as an ImageJ calibration (micron) in the file metadata, was added to each of directories under /inst/extdata/.
* R (>= 3.6.0) is required. The packages `tiff`, `rstudioapi` and `utils` are now in Imports.
* The vignette is rewritten for the current functions and renamed `Introduction_to_AiES`.
* The help pages of all functions are updated.
