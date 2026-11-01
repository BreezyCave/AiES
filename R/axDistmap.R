# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#' @title Create distance map and binary image from image file
#' @description
#' `axDistmap` processes image files (TIFF) to create distance maps and optionally binary images.
#' It also computes various image features and exports them as a text file.
#'
#' Images are resampled to a fixed physical scale (micrometers per pixel) rather than to a fixed
#' pixel width, so that data acquired on different microscopes and cameras are directly comparable.
#' Grayscale and RGB color phase-contrast images are supported via the `colorMode` argument.
#'
#' @param pixelSize Numeric or NULL. Physical size of one pixel of the *input* image, in micrometers.
#' If NULL, the file metadata (OME-XML, TIFF resolution tags) is consulted; if those are absent, scale matching
#' is disabled and `resizeW` is used instead. Metadata formats differ between manufacturers, so
#' specifying `pixelSize` explicitly is recommended. If both are available, the specified value is
#' used, and a warning reports both values. (default: NULL)
#' @param targetPixelSize Numeric or NULL. Pixel size, in micrometers, to resample the image to.
#' The default is the scale at which the built-in models were trained (0.791746 um/px).
#' To process the image without resizing it, set this parameter to the same value as pixelSize.
#' Set to NULL to ignore physical scale entirely and fall back to `resizeW`. (default: 0.791746)
#' @param resizeW Numeric or NULL. Target width in pixels, used only when the physical scale cannot
#' be resolved or when `targetPixelSize = NULL`. If NULL, the image is processed without resizing.
#' (default: 900)
#' @param subBack Numeric. Connected components (objects) whose area is at most this many pixels at the
#' reference scale (0.791746 um/px; 30 px = 18.8 um^2) are removed. With autoParam = TRUE, the value is
#' converted automatically so that the same physical area is removed at any image resolution.
#' (default: 30)
#' @param autoParam Logical. If TRUE, the pixel-unit segmentation parameters (`subBack` and the
#' local threshold window) are rescaled so that they cover a constant physical size. The Haralick
#' GLCM scales are always kept at 1 and 2 px. Has no effect when the physical scale is unknown.
#' Set to FALSE only if `targetPixelSize` matches the reference scale, or if you need to reproduce the
#' behavior of versions before 1.0.0 (together with `targetPixelSize = NULL`).
#' (default: TRUE)
#' @param binaryImg Logical. If TRUE, exports binary image files. (default: FALSE)
#' @param featureSet Character. If "standard", exports the minimum set of features required by the
#' built-in models of this package, plus Area. If "all", exports all computed features. (default: "standard").
#' @param colorMode Character. How the input image is reduced to a single grayscale channel:
#' \itemize{
#'   \item "gray": single-channel input, or channels that all contain the same grayscale image
#'         (e.g. a grayscale image stored in an RGB container). Channel 1 is used.
#'   \item "rgb": a true RGB color image, e.g. a phase-contrast image recorded with a color
#'         camera. The three channels are averaged.
#' }
#' (default: "gray")
#' @param imgType Character. Output image format: "png", "jpg", or "tiff" (default: "tiff").
#' @param folder_paths Character vector. Folder paths to process. If NULL, a folder selection dialog will be shown. (default: NULL).
#' @param output_path Character. Output folder. If NULL, a folder selection dialog will be shown. (default: tempdir()).
#'
#' @return This function doesn't return a value directly, but produces the following outputs:
#' \itemize{
#'   \item A distance map image file (format specified by `imgType`)
#'   \item A text file containing computed image features (named `<image>_<date>_ImageData.txt`)
#'   \item A parameter file (named `<folder>_<date>_Parameters.csv`)
#'   \item (Optional) A binary image file if `binaryImg = TRUE`
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Reads the input image. Single-Channel (grayscale) images and 3-channel RGB images are supported
#' 2. Reduces the image to a single grayscale channel according to `colorMode`
#' 3. Resamples the image to `targetPixelSize` (or to `resizeW` if the scale is unknown or `targetPixelSize = NULL`)
#' 4. Rescales the pixel-unit parameters if `autoParam = TRUE`
#' 5. Applies contrast adjustment, filtering and local thresholding
#' 6. Computes a distance map, then shape, Haralick, and moment features
#' 7. Exports the processed images and the computed features
#'
#' @section Input Images:
#' Axons are segmented as structures darker than their local surroundings, as in phase-contrast
#' images. Only TIFF files are processed. Images must occupy the full intensity range of their file
#' format. Data from a camera with fewer bits than the TIFF container (e.g. 12-bit data in a 16-bit
#' TIFF) should be exported as OME-TIFF, which records the significant bits (the actual bit depth of
#' the data) and is rescaled automatically, or the data should be rescaled to use the full available
#' range (e.g., 0-65535 for 16 bit images) within the acquisition software. A warning is issued
#' for suspiciously dark images. The pixel size should be known (from OME-XML or TIFF metadata, or via
#' `pixelSize`) to match images to the scale of the built-in models. If it cannot be determined,
#' the image is resized to resizeW, and a warning is issued.
#' Multi-frame TIFF files, lossy compression, contrast-adjusted images and burned-in scale bars should
#' be avoided.
#' See \code{\link{AiES-package}} (\code{?AiES}) for recommended and not recommended formats.
#'
#' @section Physical Scale:
#' Every pixel-unit parameter of the pipeline defines a physical size once the scale is fixed.
#' At the reference scale of 0.791746 um/px, `subBack = 30` corresponds to 18.8 um^2, the local
#' threshold window spans 13.5 x 13.5 um, and the Haralick GLCM offsets span 0.79 and 1.58 um.
#' Resampling to the reference scale restores all of them simultaneously, which is what
#' allows a classifier trained on one system to be applied to another.
#'
#' If the `targetPixelSize` is set to a scale different from the reference scale, parameter rescaling keeps
#' the physical size of the removed objects and of the local threshold window unchanged. However, the resulting
#' feature values will not be numerically equivalent to those obtained at the reference scale. For
#' example, since Haralick features are computed at fixed pixel scales, they will sample texture at a
#' finer physical distance. Therefore, instead of using a built-in model trained at the reference scale,
#' it is recommended to retrain the classifier using \code{\link{axSvm}}.
#'
#' @section Color Handling:
#' Depending on how a TIFF file was exported, the same field of view may be stored as a 2D grayscale
#' array, a 3-channel RGB image, or a multi-frame stack. `axDistmap` reads every channel or frame,
#' reports the detected structure, and then reduces the image to a single channel according to `colorMode`.
#' "gray" (default) uses the first channel, while "rgb" averages the three channels. If "rgb" is selected
#' but the image does not have exactly three channels, the file is skipped with a warning. Multi-frame
#' files (z-stacks, time series) are not recommended, as only the first frame is analyzed.
#'
#' @section Feature Export:
#' All features are exported under short names (see section "Feature names"), which are the names used
#' by the built-in models. If featureSet = "standard", the predictors of the built-in models (Per, Ras, Ran,
#' ecc, mjx, ENT, SVA, SV2, CON and cir) and Area (the area weight for AII/DI) are exported. The 'cir'
#' feature represents Circularity and is calculated as:
#' cir = (Area * pi * 4) / (Per^2)
#'
#' If featureSet = "all", all computed features from EBImage's computeFeatures functions
#' (shape, moment, and haralick) plus 'cir' are exported.
#'
#' @section Parameter file:
#' For every output folder a file `<folder>_<date>_Parameters.csv` is written next to the
#' feature files, with one row per processed image: `FileName` (path relative to the input
#' folder), `DataFile` (the matching ImageData.txt), `Date`, `AiESVersion`, `WidthPx`,
#' `HeightPx`, `Channels`, `ContainerBits`, `SignificantBits`, `BitRescaled`, `ColorMode`,
#' `ArgPixelSize` (the `pixelSize` argument), `MetaPixelSize` and `MetaSource` (value found in
#' the file and where), `PixelSize` (effective scale, um/px, of the analyzed image; NA when
#' unknown), `ScaleMode` (`scale-matched`, `fixed-width` or `no-resize`), `ResizeW` (working
#' width, px), `SubBack` and `ThreshWH` (values actually used) and `FeatureSet`.
#' \code{\link{axSvm}} reads `PixelSize` from this file to record the training scale with the
#' model. The feature files themselves contain features only.
#'
#' @section Feature names:
#' Each entry reads: short name = original EBImage column name (cir is calculated by AiES).
#' \itemize{
#'   \item Shape: Area = s.area (area weight for AII / DI), Per = s.perimeter,
#'         Ram = s.radius.mean, Ras = s.radius.sd, Ran = s.radius.min, Rax = s.radius.max.
#'   \item Moment: mcx = m.cx, mcy = m.cy, mjx = m.majoraxis, ecc = m.eccentricity, the = m.theta.
#'   \item Haralick, scale 1: ASM = h.asm.s1, CON = h.con.s1, COR = h.cor.s1, VAR = h.var.s1,
#'         IDM = h.idm.s1, SAV = h.sav.s1, SVA = h.sva.s1, SEN = h.sen.s1, ENT = h.ent.s1,
#'         DVA = h.dva.s1, DEN = h.den.s1, F12 = h.f12.s1, F13 = h.f13.s1.
#'   \item Haralick, scale 2: SM2 = h.asm.s2, ON2 = h.con.s2, OR2 = h.cor.s2, AR2 = h.var.s2,
#'         DM2 = h.idm.s2, AV2 = h.sav.s2, SV2 = h.sva.s2, SE2 = h.sen.s2, NT2 = h.ent.s2,
#'         VA2 = h.dva.s2, EN2 = h.den.s2, F22 = h.f12.s2, F32 = h.f13.s2.
#'   \item Circularity: cir.
#' }
#'
#' @note
#' - Output files are saved under `output_path` in a dated subfolder.
#' - Feature columns are exported under the short names above. Text files written by earlier
#'   versions (original EBImage names) are still accepted by \code{\link{axSvm}} and
#'   \code{\link{axQnt}} and converted automatically.
#' - Haralick features are computed at scales of 1 and 2 px (scale 1: ASM, CON, ...; scale 2: SM2, ON2, ...)
#'   at any resolution.
#'
#' @examples
#'
#' # Basic usage with default parameters
#' # NOTE: This example requires a GUI environment for interactive folder selection.
#' if (interactive()){
#' axDistmap()
#' }
#'
#' \donttest{
#' # Declare the scale of the input images explicitly
#' # (the sample images are already resized to 900 px width, i.e. 0.791746 um/px)
#' img_dir <- system.file("extdata", "Degenerate_Images", package = "AiES")
#' axDistmap(pixelSize = 0.791746, folder_paths = img_dir, output_path = tempdir())
#'
#' # Work at the input's own scale (parameters are rescaled automatically; for the
#' # sample images this equals the reference scale)
#' axDistmap(pixelSize = 0.791746, targetPixelSize = 0.791746,
#' folder_paths = img_dir, output_path = tempdir())
#'
#' # Average the three channels (the sample images store a grayscale image in an RGB container)
#' axDistmap(pixelSize = 0.791746, colorMode = "rgb",
#' folder_paths = img_dir, output_path = tempdir())
#'
#' # Legacy behavior: ignore physical scale, resize to a fixed width
#' axDistmap(targetPixelSize = NULL, resizeW = 900, featureSet = "all",
#' folder_paths = img_dir, output_path = tempdir())
#' }
#'
#' # Process with a custom width and background threshold
#' img_dir <- system.file("extdata", "Degenerate_Images", package = "AiES")
#' axDistmap(subBack = 20, targetPixelSize = NULL, resizeW = 300, imgType = "png",
#' folder_paths = img_dir, output_path = tempdir())
#'
#'
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage Image
#' @importFrom fs dir_ls
#' @importFrom utils write.table
#' @importFrom stats na.omit setNames
#'
#' @export
axDistmap <- function(pixelSize = NULL,
                      targetPixelSize = 0.791746,
                      resizeW = 900,
                      subBack = 30,
                      autoParam = TRUE,
                      binaryImg = FALSE,
                      featureSet = "standard",
                      colorMode = "gray",
                      imgType = "tiff",
                      folder_paths = NULL,
                      output_path = tempdir()){

    ###imgType check
    if(!(imgType %in% c("tiff","png","jpg"))) imgType <- "tiff"
    ## EBImage::writeImage() accepts "jpeg", not "jpg"; keep ".jpg" as the extension
    imgTypeOut <- if (imgType == "jpg") "jpeg" else imgType

    ###featureSet / colorMode check
    featureSet <- match.arg(featureSet, c("standard", "all"))
    colorMode  <- match.arg(colorMode, c("gray", "rgb"))

    sdate <- Sys.Date()
    ###############Definition of function######################
    ##Function 1
    ## Functions are executed depending on the user environment
    select_folder <- function(caption = "Select folder") {
        if (rstudioapi::isAvailable()) {
            return(rstudioapi::selectDirectory(caption = caption))
        } else if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.dir(caption = caption))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_dir(title = caption)$res)
        } else {
            message("No GUI folder selection method available")
            return(NULL)
        }
    }

    # Selection or specification of folders
    if (is.null(folder_paths)) {
        folder_paths <- character(0)
        repeat {
            selected_folder <- select_folder()
            if (is.null(selected_folder)) {
                message("Folder selection cancelled.")
                break
            }
            if (!(selected_folder %in% folder_paths)) {
                folder_paths <- c(folder_paths, selected_folder)
                message(sprintf("Added folder: %s", selected_folder))
            } else {
                message(sprintf("Folder already selected: %s. Skipping.", selected_folder))
            }
            if (tolower(readline(prompt = "Select another folder? (y/n): ")) != "y") break
        }
        if (length(folder_paths) == 0) {
            message("No folders selected. Exiting function.")
            return(NULL)
        }
    }

    folder_paths <- as.list(folder_paths)
    folder_paths <- Filter(dir.exists, folder_paths)
    folder_paths <- unique(folder_paths)

    if (length(folder_paths) == 0) {
        message("No valid folders to process. Exiting function.")
        return(NULL)
    }

    # Create output folder
    if (is.null(output_path)) {
        output_path <- select_folder("Select output folder")
        if (is.null(output_path) || output_path == "") output_path <- getwd()
    }
    if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


    # Process each folder
    results <- lapply(folder_paths, function(folder_path) {
        if (!dir.exists(folder_path)) {
            message(sprintf("Specified folder does not exist: %s. Skipping.", folder_path))
            return(NULL)
        }

        tiff_files <- fs::dir_ls(folder_path, regexp = "\\.tiff?$", recurse = TRUE)

        if (length(tiff_files) == 0) {
            message(sprintf("No TIFF files found in the folder: %s. Skipping.", folder_path))
            return(NULL)
        }

        process_file <- function(file) {
            message(sprintf("Processing file: %s", file))

            output_dir <- file.path(output_path,
                                    paste0(basename(folder_path), "_", sdate, "_output_files"))
            if (!dir.exists(output_dir)) dir.create(output_dir)

            ## Read, reduce to one channel, resample (shared implementation)
            grayImage <- .aiesReadGray(file, colorMode,
                                       pixelSize, targetPixelSize, resizeW)
            if (is.null(grayImage)) return(NULL)
            scale <- attr(grayImage, "aiesScale")
            info  <- attr(grayImage, "aiesInfo")

            ## Rescale the pixel-unit parameters to the working resolution
            prm <- .aiesScaleParams(scale$k, subBack, autoParam)
            if (!isTRUE(all.equal(scale$k, 1))) {
                message(sprintf("  Parameters (k = %.3f): subBack = %s, thresh = %s",
                                scale$k, prm$subBack, prm$threshWH))
            }

            ## Segment and compute features (shared implementation)
            seg <- .aiesSegment(grayImage, prm$subBack, prm$threshWH, prm$haralickScales)
            singleData <- seg$features
            bnrySeg <- seg$segmentation
            nDm <- seg$distmap
            invisible({rm(grayImage, seg); gc()})

            if (featureSet == "standard") {
                singleData <- singleData[, .aiesStandardColumns, drop = FALSE]
            }
            dataFile <- paste0(basename(file), "_", sdate, "_ImageData.txt")
            prmRow <- .aiesParamRow(.aiesRelPath(file, folder_path), dataFile, sdate,
                                    info, colorMode, scale, prm, featureSet)

            #####export data file
            write.table(singleData,
                        file.path(output_dir, dataFile),
                        sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
            message(dataFile)

            if (binaryImg == TRUE) {
                dmrmabw <- 1*(bnrySeg > 0)
                file.name <- file.path(output_dir, paste0(basename(file), "_", sdate, "_Binary.", imgType))
                writeImage(dmrmabw, file.name, type = imgTypeOut, quality = 100)
                message(file.name)
            }

            nDmrmadm <- nDm*(bnrySeg > 0)   # distance map

            file.name <- file.path(output_dir, paste0(basename(file), "_", sdate, "_DistMap.", imgType))
            writeImage(nDmrmadm, file.name, type = imgTypeOut, quality = 100)
            message(file.name)
            prmRow
        }

        folder_results <- lapply(tiff_files, process_file)

        ## Parameter record: one row per processed image
        prmRows <- do.call(rbind, unname(folder_results))
        if (!is.null(prmRows) && nrow(prmRows) > 0) {
            output_dir <- file.path(output_path,
                                    paste0(basename(folder_path), "_", sdate, "_output_files"))
            prmFile <- file.path(output_dir,
                                 paste0(basename(folder_path), "_", sdate, "_Parameters.csv"))
            utils::write.csv(prmRows, prmFile, row.names = FALSE, na = "")
            message(prmFile)
        }

        message(sprintf("All files in folder %s processed successfully.", folder_path))

        return(folder_results)
    })

    message("All folders processed successfully.")
}
