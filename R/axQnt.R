# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#' @title Quantify Axon Integrity and Degeneration Indices from Images or Feature Data
#' @description
#' `axQnt` calculates the Axon Integrity Index (AII) and Degeneration Index (DI) from axon image files (.tiff) or from text files (.txt) of image features exported by \code{\link{axDistmap}}.
#' The function uses a pre-trained SVM model to classify axonal regions and computes indices based on classified areas.
#' `axQnt` supports both interactive GUI-based file/folder selection and direct path specification for input and output, with cross-platform compatibility.
#'
#' Either one of the built-in models
#' (`model = "g0.1_APA"` or `"g0.2_ACC"`) or any model produced by \code{\link{axSvm}}
#' can be used.
#'
#' @param model Character or NULL. Built-in model to apply: "g0.1_APA" or "g0.2_ACC".
#'        Ignored when `svm_model_path` is given. If both `model` and `svm_model_path` are NULL,
#'        a file selection dialog is shown to choose one model file. (default: "g0.1_APA")
#' @param imprtImg Logical. If TRUE, imports and processes image files (.tiff). If FALSE, uses precomputed feature data (.txt). (default: TRUE)
#'        The image-processing arguments (`subBack`, `pixelSize`, `targetPixelSize`, `resizeW`, `autoParam`,
#'        `colorMode`, `binaryImg` and `distMapImg`) are used only when `imprtImg = TRUE`. Arguments shared
#'        with \code{\link{axDistmap}} are described identically in both functions.
#' @param subBack Numeric. Connected components (objects) whose area is at most this many pixels at the
#' reference scale (0.791746 um/px; 30 px = 18.8 um^2) are removed. With autoParam = TRUE, the value is
#' converted automatically so that the same physical area is removed at any image resolution.
#' (default: 30)
#' @param pixelSize Numeric or NULL. Physical size of one pixel of the *input* image, in micrometers.
#' If NULL, the file metadata (OME-XML, TIFF resolution tags) is consulted; if those are absent, scale matching
#' is disabled and `resizeW` is used instead. Metadata formats differ between manufacturers, so
#' specifying `pixelSize` explicitly is recommended. If both are available, the specified value is
#' used, and a warning reports both values. (default: NULL)
#' @param targetPixelSize Numeric or NULL. Pixel size, in micrometers, to resample the image to.
#' The default is the scale recorded with the model (0.791746 um/px for the built-in models, and for
#' models without a recorded scale), so that images are brought to the scale the model was trained at.
#' To process the image without resizing it, set this parameter to the same value as pixelSize.
#' Set to NULL to ignore physical scale entirely and fall back to `resizeW`.
#' (default: taken from the model; leave the argument out to use it)
#' @param resizeW Numeric or NULL. Target width in pixels, used only when the physical scale cannot
#' be resolved or when `targetPixelSize = NULL`. If NULL, the image is processed without resizing.
#' (default: 900)
#' @param autoParam Logical. If TRUE, the pixel-unit segmentation parameters (`subBack` and the
#' local threshold window) are rescaled so that they cover a constant physical size. The Haralick
#' GLCM scales are always kept at 1 and 2 px. Has no effect when the physical scale is unknown.
#' Set to FALSE only if `targetPixelSize` matches the reference scale, or if you need to reproduce the
#' behavior of versions before 1.0.0 (together with `targetPixelSize = NULL`).
#' (default: TRUE)
#' @param colorMode Character. How the input image is reduced to a single grayscale channel:
#' \itemize{
#'   \item "gray": single-channel input, or channels that all contain the same grayscale image
#'         (e.g. a grayscale image stored in an RGB container). Channel 1 is used.
#'   \item "rgb": a true RGB color image, e.g. a phase-contrast image recorded with a color
#'         camera. The three channels are averaged.
#' }
#' A single-channel image is processed as is, whatever colorMode is. (default: "gray")
#' @param binaryImg Logical. If TRUE, exports binary image files. (default: FALSE)
#' In `axQnt`, the files are written as TIFF.
#' @param distMapImg Logical. If TRUE, exports the distance map image of each file (the same image as
#'        written by \code{\link{axDistmap}}), as a TIFF file. Used with `imprtImg = TRUE`. (default: FALSE)
#' @param expSip Logical. If TRUE, exports single image prediction results as .csv files. (default: TRUE)
#' @param svm_model_path Character or NULL. Path to one user-supplied SVM model file (.RData or .svm),
#'        typically produced by \code{\link{axSvm}} (not a folder). The file may contain other objects (e.g.
#'        a saved workspace); the first object of class "svm" is used, under any name, and the other objects
#'        do not affect the analysis. If given, it takes precedence over `model`. (default: NULL)
#' @param input_dirs Character. Path to the folders containing input files (.tiff or .txt). If NULL, a folder selection dialog will be shown. (default: NULL)
#' @param output_dir Character. Path to the output directory. If NULL, a folder selection dialog will be shown; if empty, current working directory is used. One subfolder per input folder is created below it (see section "Output Folders"). (default: tempdir())
#'
#' @return Invisibly returns the combined summary table (NULL if no file could be processed) and
#' generates the following outputs:
#' \itemize{
#'   \item `Summary_<date>.csv` in `output_dir`: the calculated Axon Integrity Index (AII) and Degeneration Index (DI) for every file of every input folder. `FileName` is the path relative to the input folder. Not written if no file could be processed.
#'   \item In the subfolder of each input folder: `<folder>_<date>_Summary.csv` (the same table for that folder only).
#'   \item (Optional) Single image prediction results as .csv files if `expSip = TRUE`.
#'   \item (Optional) Binary images if `binaryImg = TRUE`, and distance map images (`<file>_<date>_DistMap.tiff`) if `distMapImg = TRUE`.
#'   \item (Image route only) In the same subfolder, `<folder>_<date>_Parameters.csv`: one row per processed image, with the same columns as the parameter file of \code{\link{axDistmap}} (see there). `FileName` is the path relative to the input folder; `DataFile` is the single image prediction file, or empty if `expSip = FALSE`; `FeatureSet` is "all".
#' }
#'
#' @details
#' The function implements a complete quantification pipeline:
#' \enumerate{
#'   \item \strong{Model loading}: Loads a built-in model, or a pre-trained SVM model from the specified file or via GUI selection. Automatically detects SVM objects within the loaded file.
#'   \item \strong{Predictor detection}: Reads the predictors the model was trained on directly from the model object, and maps them back to the corresponding \code{\link{axDistmap}} column names. Any model produced by \code{\link{axSvm}} therefore works without further configuration.
#'   \item \strong{Input handling}: If input paths are not specified, the user is prompted to select folders via a GUI dialog.
#'   \item \strong{Data preprocessing}: Depending on \code{imprtImg}, either processes image files (.tiff) using exactly the same pipeline as \code{\link{axDistmap}} or directly utilizes precomputed feature data (.txt).
#'   \item \strong{SVM classification}: Classifies axonal regions into "Degenerate" or "Intact".
#'   \item \strong{Index calculation}: Calculates Axon Integrity Index (AII) and Degeneration Index (DI) based on classified areas:
#'         \deqn{AII = \frac{\text{Area of Intact Axons}}{\text{Total Axonal Area}}}
#'         \deqn{DI = \frac{\text{Area of Degenerate Axons}}{\text{Total Axonal Area}}}
#'   \item \strong{Result saving}: Exports summary and optional detailed prediction files with unique filenames to prevent overwriting.
#' }
#'
#' @section Input Images:
#' Axons are segmented as structures darker than their local surroundings, as in phase-contrast
#' images. Only TIFF files are processed. Images must occupy the full intensity range of their file
#' format. Data from a camera with fewer bits than the TIFF container (e.g. 12-bit data in a 16-bit
#' TIFF) should be exported as OME-TIFF, which records the significant bits and is rescaled
#' automatically, or be scaled to the full range in the acquisition software. A warning is issued
#' for suspiciously dark images. The pixel size must be known (from OME-XML or TIFF metadata, or via
#' `pixelSize`) for the built-in models to apply.
#' Lossy compression, contrast-adjusted images and burned-in scale bars should be avoided.
#' See \code{\link{AiES-package}} (\code{?AiES}) for recommended and not recommended formats.
#'
#' @section Physical Scale:
#' The built-in models were trained on images analyzed at 0.791746 um/px. By default, `axQnt`
#' resamples every input image to match the model's scale. This keeps pixel-based pipeline
#' parameters such as background area threshold, local threshold window, and Haralick GLCM
#' offsets consistent with the physical distances used during training.
#'
#' If an image's scale differs from that of the model (for example, when pixel size is unknown
#' and the image is resized with `resizeW`, or when `targetPixelSize` is set), a warning at the
#' end summarizes the number of such images, their working scales, and the ratio to the model's
#' scale. The working scale for each image is also saved in the `PixelSize` column.
#'
#' Features such as Perimeter, Circularity, and the Haralick texture features are sensitive to image
#' resolution. Changes in pixel size affect these values, and simply rescaling does not fully
#' correct for this. For accurate results, create a new model with \code{\link{axSvm}} using data acquired
#' at the new resolution.
#'
#' @section Output Folders:
#' For each input folder, a subfolder `<folder>_<date>_output_files` is created in `output_dir`
#' (the same naming as \code{\link{axDistmap}}). If several input folders have the same name, they are
#' numbered (`<folder>_2_<date>_output_files`, ...). The single image prediction files, binary images,
#' the folder's summary and its parameter file are written there. The combined `Summary_<date>.csv`
#' of all input folders is written directly to `output_dir`.
#'
#' @section File Naming and Overwrite Policy:
#' \itemize{
#'   \item Output files are automatically named with timestamps in ISO 8601 format (YYYY-MM-DD).
#'   \item If a file name already exists, a unique name with a numeric suffix is generated to avoid overwriting.
#'   \item If the output directory does not exist, it will be created automatically.
#' }
#'
#' @section GUI Support:
#' \itemize{
#'   \item RStudio (rstudioapi), tcltk, and svDialogs are supported for file/folder selection.
#'   \item If no GUI is available, the user is prompted to enter the path manually.
#'   \item Cross-platform compatibility for Windows, macOS, and Linux.
#' }
#'
#' @note
#' \itemize{
#'   \item Input .txt files are expected to carry the short feature names exported by
#'         \code{\link{axDistmap}}. Files written by earlier versions (original EBImage names)
#'         are converted automatically.
#'   \item The column Area must be present, as it provides the area weighting for AII and DI.
#'   \item Output files include timestamp in ISO 8601 format (YYYY-MM-DD).
#' }
#'
#' @examples
#' # Interactive mode: process .tiff images with GUI dialogs
#' # NOTE: This example requires a GUI environment for interactive folder selection.
#' if (interactive()){
#' axQnt(imprtImg = TRUE, expSip = TRUE)
#' }
#'
#' \donttest{
#' # Apply the built-in g0.1 APA model to package-included images
#' img_dir <- system.file("extdata", "Degenerate_Images", package = "AiES")
#' axQnt(model = "g0.1_APA", imprtImg = TRUE, input_dirs = img_dir, output_dir = tempdir())
#'
#' # Declare the scale of the input images explicitly
#' axQnt(model = "g0.1_APA", pixelSize = 0.791746, imprtImg = TRUE,
#'       input_dirs = img_dir, output_dir = tempdir())
#'
#' # Apply the built-in g0.2 ACC model instead
#' axQnt(model = "g0.2_ACC", imprtImg = TRUE, input_dirs = img_dir, output_dir = tempdir())
#'
#' # Apply a model re-trained with axSvm()
#' my_model <- file.path(tempdir(), "my_classifier.svm")
#' if (file.exists(my_model)) {
#'   axQnt(svm_model_path = my_model, imprtImg = TRUE,
#'         input_dirs = img_dir, output_dir = tempdir())
#' }
#' }
#'
#' # Utilize package-included precomputed feature data from .txt files
#' txt_dir <- system.file("extdata", "Degenerate_txt", package = "AiES")
#' axQnt(model = "g0.1_APA", imprtImg = FALSE, input_dirs = txt_dir, output_dir = tempdir())
#'
#'
#' @importFrom data.table fread fwrite rbindlist data.table
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage Image
#' @importFrom stats predict na.omit setNames
#'
#' @export
axQnt <- function(model = "g0.1_APA",
                  imprtImg = TRUE,
                  subBack = 30,
                  pixelSize = NULL,
                  targetPixelSize = NULL,
                  resizeW = 900,
                  autoParam = TRUE,
                  colorMode = "gray",
                  binaryImg = FALSE,
                  distMapImg = FALSE,
                  expSip = TRUE,
                  svm_model_path = NULL,
                  input_dirs = NULL,
                  output_dir = tempdir()){

    sdate <- Sys.Date()

    colorMode <- match.arg(colorMode, c("gray", "rgb"))

    ###############Definition of function######################
    ##Function 1
    select_folder <- function(caption = "Select folder") {
        if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
            return(rstudioapi::selectDirectory(caption = caption))
        } else if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.dir(caption = caption))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_dir(title = caption)$res)
        } else {
            cat("Enter folder path: ")
            return(readline())
        }
    }

    ##Function 2 (one file only; used to choose the model file)
    select_file <- function(caption = "Select a file") {
        if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.files(caption = caption, multi = FALSE))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_open(title = caption, multiple = FALSE)$res)
        } else {
            cat("Enter file path: ")
            return(readline())
        }
    }

    ##Function 3
    generate_unique_filename <- function(filepath) {
        if (!file.exists(filepath)) return(filepath)
        base <- sub("(\\.[^.]*)$", "", filepath)
        ext <- sub("^.*(\\.[^.]*)$", "\\1", filepath)
        i <- 1
        new_filepath <- paste0(base, "_", i, ext)
        while (file.exists(new_filepath)) {
            i <- i + 1
            new_filepath <- paste0(base, "_", i, ext)
        }
        return(new_filepath)
    }

    # 1. Resolve and load the SVM model
    if (is.null(svm_model_path)) {
        if (!is.null(model)) {
            model <- match.arg(model, names(.aiesModels))
            svm_model_path <- system.file("extdata", .aiesModels[[model]]$file, package = "AiES")
            if (svm_model_path == "") {
                stop("Built-in model file not found for '", model, "'.")
            }
        } else {
            svm_model_path <- select_file("Select an SVM model file")
            if (length(svm_model_path) != 1 || is.na(svm_model_path) || svm_model_path == "") {
                return(message("SVM model not selected."))
            }
        }
    }

    ## The path must be one existing file (checked after resolution, since
    ## svm_model_path is NULL when a built-in model or the dialog is used)
    if (length(svm_model_path) != 1 || is.na(svm_model_path)) {
        stop("`svm_model_path` must be the path to one model file.", call. = FALSE)
    }
    if (dir.exists(svm_model_path)) {
        stop("`svm_model_path` must be a model file (.svm or .RData), not a folder:\n", svm_model_path, call. = FALSE)
    }
    if (!file.exists(svm_model_path)) {
        stop("Model file not found:\n", svm_model_path, call. = FALSE)
    }

    ## Load into a separate environment, so that objects saved in the file
    ## (e.g. a whole workspace) cannot overwrite the arguments of this function
    model_env <- new.env()
    loaded_objects <- load(svm_model_path, envir = model_env)
    svm_model <- NULL
    for (obj_name in loaded_objects) {
        obj <- get(obj_name, envir = model_env, inherits = FALSE)
        if (inherits(obj, "svm")) {
            svm_model <- obj
            break
        }
    }
    rm(model_env)   # release the other objects of the file
    if (is.null(svm_model)) {
        stop("No SVM model found in the loaded file:\n", svm_model_path)
    }

    ## Predictors the model was actually trained on (short names)
    ftrSvm <- .aiesModelFeatures(svm_model)
    if (is.null(ftrSvm) || length(ftrSvm) == 0) {
        stop("Could not determine the predictors of the loaded SVM model.")
    }
    model_label <- attr(svm_model, "aiesModel")
    message("Model: ", if (is.null(model_label)) basename(svm_model_path) else model_label,
            " | features: ", paste(ftrSvm, collapse = ", "))

    ## Scale the model was trained at: explicit argument > model attribute > preset > reference
    modelPixelSize <- attr(svm_model, "aiesPixelSize")
    if (is.null(modelPixelSize) && !is.null(model_label) && model_label %in% names(.aiesModels)) {
        modelPixelSize <- .aiesModels[[model_label]]$targetPixelSize
    }
    if (missing(targetPixelSize)) {
        targetPixelSize <- if (is.null(modelPixelSize)) .aiesRefPixelSize else modelPixelSize
    }
    if (!is.null(targetPixelSize)) {
        message("  target scale: ", sprintf("%.6f", targetPixelSize), " um/px")
    }

    # 2. Select input directory
    if (is.null(input_dirs)) {
        input_dirs <- character(0)
        continue_selection <- "y"
        repeat {
            selected_folder <- select_folder(if (imprtImg) "Select folder containing TIFF images" else "Select folder containing TXT feature files")
            if (is.null(selected_folder)) {
                message("Folder selection cancelled.")
                break
            }
            if (!(selected_folder %in% input_dirs)) {
                input_dirs <- c(input_dirs, selected_folder)
                message(sprintf("Added folder: %s", selected_folder))
            } else {
                message(sprintf("Folder already selected: %s. Skipping.", selected_folder))
            }
            continue_selection <- readline(prompt = "Select another folder? (y/n): ")
            if (tolower(continue_selection) != "y") break
        }
        if (length(input_dirs) == 0) {
            message("No input folder selected. Exiting function.")
            return(NULL)
        }
    }

    # 3. Create output folder
    if (is.null(output_dir)) {
        output_dir <- select_folder("Select output folder")
        if (is.null(output_dir) || output_dir == "") output_dir <- getwd()
    }
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    ## One output subfolder per input folder (same naming as axDistmap();
    ## input folders with the same name are numbered)
    dir_base  <- basename(sub("[/\\\\]+$", "", input_dirs))
    dir_n     <- stats::ave(seq_along(dir_base), dir_base, FUN = seq_along)
    dir_label <- ifelse(dir_n == 1L, dir_base, paste0(dir_base, "_", dir_n))
    out_dirs  <- file.path(output_dir, paste0(dir_label, "_", sdate, "_output_files"))

    # 4. Get file list (multiple folders and subfolders supported)
    file_pattern <- if (imprtImg) "\\.tiff?$" else "\\.txt$"
    file_list <- lapply(input_dirs, function(dir) {
        list.files(dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
    })
    input_files <- unlist(file_list)
    file_dir    <- rep(seq_along(input_dirs), lengths(file_list))   # input folder of each file
    if (length(input_files) == 0) return(message("No input files found."))
    for (d in unique(file_dir)) dir.create(out_dirs[d], recursive = TRUE, showWarnings = FALSE)


    # 5. Process images or feature tables
    results_list <- lapply(seq_along(input_files), function(i) {

        file <- input_files[i]
        od   <- out_dirs[file_dir[i]]      # output subfolder of this file's input folder
        info <- NULL
        if (imprtImg) {
            message(sprintf("Processing file: %s", file))

            ## Identical pipeline to axDistmap() (shared implementation)
            grayImage <- .aiesReadGray(file, colorMode,
                                       pixelSize, targetPixelSize, resizeW)
            if (is.null(grayImage)) return(NULL)
            scale <- attr(grayImage, "aiesScale")
            info  <- attr(grayImage, "aiesInfo")

            ## The working scale is compared with the model's scale after all files
            ## are processed (one summary warning; see below).

            prm <- .aiesScaleParams(scale$k, subBack, autoParam)
            seg <- .aiesSegment(grayImage, prm$subBack, prm$threshWH, prm$haralickScales)
            singleData <- seg$features

            if (binaryImg) {
                dmrmabw <- 1 * (seg$segmentation > 0)
                bin_outfile <- file.path(od, paste0(basename(file), "_", sdate, "_Binary.tiff"))
                bin_outfile <- generate_unique_filename(bin_outfile)
                EBImage::writeImage(dmrmabw, bin_outfile, type = "tiff")
            }
            if (distMapImg) {
                dm_outfile <- file.path(od, paste0(basename(file), "_", sdate, "_DistMap.tiff"))
                dm_outfile <- generate_unique_filename(dm_outfile)
                EBImage::writeImage(seg$distmap * (seg$segmentation > 0), dm_outfile, type = "tiff")
            }
            invisible({rm(grayImage, seg); gc()})

        } else {
            singleData <- tryCatch(data.table::fread(file), error = function(e) NULL)
            if (is.null(singleData)) return(NULL)
            singleData <- as.data.frame(singleData)
        }

        ## Short names (converts only files written by earlier versions, which carry EBImage names)
        singleData <- .aiesToShort(singleData)

        ## Availability checks
        ## The data carry short names; a model trained with original EBImage names
        ## (e.g. one made by v0.99.6) is matched through the same name mapping.
        ftrData <- .aiesAsShort(ftrSvm)
        missing_ftr <- setdiff(ftrData, colnames(singleData))
        if (length(missing_ftr) > 0) {
            message("Required features missing in file: ", file,
                    "\n  Missing: ", paste(missing_ftr, collapse = ", "))
            return(NULL)
        }
        if (!("Area" %in% colnames(singleData))) {
            message("Column 'Area' (area weight for AII/DI) missing in file: ", file)
            return(NULL)
        }

        # SVM classification
        newData <- singleData[, ftrData, drop = FALSE]
        colnames(newData) <- ftrSvm          # names the model was trained with
        pred <- predict(svm_model, newData, type = "class")

        # Calculate AII/DI
        intact_area <- sum((pred == "Intact") * singleData$Area)
        deg_area    <- sum((pred == "Degenerate") * singleData$Area)
        total_area  <- intact_area + deg_area
        AII <- if (total_area > 0) intact_area / total_area else NA
        DI  <- if (total_area > 0) deg_area / total_area else NA

        # Save single image prediction results
        sipName <- NA_character_
        if (expSip) {
            sip_outfile <- file.path(od, paste0(basename(file), "_", sdate, "_SIP.csv"))
            sip_outfile <- generate_unique_filename(sip_outfile)
            data.table::fwrite(data.frame(Pred = pred, singleData), file = sip_outfile, sep = ",")
            sipName <- basename(sip_outfile)
        }

        ## File name relative to its input folder (same as FileName of the parameter file)
        relName <- .aiesRelPath(file, input_dirs[file_dir[i]])

        ## Parameter record and working scale (image route only)
        prmRow <- if (is.null(info)) NULL else
            .aiesParamRow(relName, sipName, sdate, info, colorMode, scale, prm, "all")
        workScale <- if (is.null(info)) NA_real_ else scale$effPixelSize

        return(list(sum = data.table::data.table(FileName = relName,
                                                 AxonIntegrityIndex = AII,
                                                 DegenerationIndex = DI),
                    par = prmRow,
                    scl = workScale))
    })

    ## Nothing processed: no summary files
    if (all(vapply(results_list, is.null, logical(1)))) {
        message("No file could be processed. No summary file was written.")
        return(invisible(NULL))
    }

    ## Working scale vs. the model's scale: one warning for all images. Differences
    ## below 0.01% are rounding only (e.g. fixed-width resizing) and are not reported.
    if (!is.null(modelPixelSize)) {
        scl <- vapply(results_list, function(r) if (is.null(r)) NA_real_ else r$scl, numeric(1))
        scl <- scl[is.finite(scl)]
        mis <- scl[abs(scl / modelPixelSize - 1) > 1e-4]
        if (length(mis) > 0) {
            tab <- table(round(mis, 4))
            u   <- as.numeric(names(tab))
            warning(sprintf(
                "%d of %d images were analyzed at a working scale different from the model's %.4f um/px: %s. See PixelSize in the parameter file.",
                length(mis), length(scl), modelPixelSize,
                paste(sprintf("%.4f um/px (x%.3f, %d images)", u, u / modelPixelSize,
                              as.integer(tab)), collapse = "; ")),
                call. = FALSE)
        }
    }

    # Per input folder: summary and parameter record, in its output subfolder
    for (k in unique(file_dir)) {
        idx   <- which(file_dir == k)
        res_k <- data.table::rbindlist(lapply(results_list[idx], `[[`, "sum"), fill = TRUE)
        if (nrow(res_k) == 0) next
        sm_k <- generate_unique_filename(
            file.path(out_dirs[k], paste0(dir_label[k], "_", sdate, "_Summary.csv")))
        data.table::fwrite(res_k, file = sm_k, sep = ",")
        message("Summary file saved: ", sm_k)

        prm_k <- data.table::rbindlist(lapply(results_list[idx], `[[`, "par"), fill = TRUE)
        if (nrow(prm_k) > 0) {
            pf_k <- generate_unique_filename(
                file.path(out_dirs[k], paste0(dir_label[k], "_", sdate, "_Parameters.csv")))
            data.table::fwrite(prm_k, file = pf_k, sep = ",", na = "")
            message("Parameter file saved: ", pf_k)
        }
    }

    # Combined summary of all input folders
    results_summary <- data.table::rbindlist(lapply(results_list, `[[`, "sum"), fill = TRUE)
    summary_outfile <- file.path(output_dir, paste0("Summary_", sdate, ".csv"))
    summary_outfile <- generate_unique_filename(summary_outfile)
    data.table::fwrite(results_summary, file = summary_outfile, sep = ",")
    message("Summary file saved (all folders): ", summary_outfile)
    invisible(results_summary)


}#####end of function
