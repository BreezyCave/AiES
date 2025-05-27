# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#' @title Quantify Axon Integrity and Degeneration Indices from Images or Feature Data
#' @description
#' `axQnt` calculates the Axon Integrity Index (AII) and Degeneration Index (DI) from axon image files (.tiff) or precomputed feature data (.txt).
#' The function uses a pre-trained SVM model to classify axonal regions and computes indices based on classified areas.
#' `axQnt` supports both interactive GUI-based file/folder selection and direct path specification for input and output, with cross-platform compatibility.
#'
#'
#' @param imprtImg Logical. If TRUE, imports and processes image files (.tiff). If FALSE, uses precomputed feature data (.txt). (default: TRUE)
#' @param subBack Numeric. Area threshold in pixels: connected components (objects) with area less than or equal to this value will be removed as background. (default: 30)
#' @param resizeW Numeric. Target width in pixels for resizing each imported image (default: 900). Larger values may improve analysis accuracy up to a point, but also increase memory and computation time. Excessively large values may not further improve results.
#' @param binaryImg Logical. If TRUE, exports binary image files. (default: FALSE)
#' @param expSip Logical. If TRUE, exports single image prediction results as .csv files. (default: TRUE)
#' @param svm_model_path Character. Path to a pre-trained SVM model file (.RData or .svm). If NULL, a file selection dialog will be shown. (default: NULL)
#' @param input_dirs Character. Path to the folders containing input files (.tiff or .txt). If NULL, a folder selection dialog will be shown. (default: NULL)
#' @param output_dir Character. Path to the output directory. If NULL, a folder selection dialog will be shown; if empty, current working directory is used. (default: NULL)
#'
#' @return The function does not return values directly but generates the following outputs:
#' \itemize{
#'   \item A summary .csv file containing the calculated Axon Integrity Index (AII) and Degeneration Index (DI) for each file.
#'   \item (Optional) Single image prediction results as .csv files if `expSip = TRUE`.
#'   \item (Optional) Processed image data as .txt files if `imprtImg = TRUE`.
#' }
#'
#' @details
#' The function implements a complete quantification pipeline:
#' \enumerate{
#'   \item \strong{Model loading}: Loads a pre-trained SVM model from the specified file or via GUI selection. Automatically detects SVM objects within the loaded file.
#'   \item \strong{Input handling}: If input paths are not specified, the user is prompted to select folders via a GUI dialog (cross-platform support for RStudio, tcltk, svDialogs, or manual input).
#'   \item \strong{Data preprocessing}: Depending on \code{imprtImg}, either processes image files (.tiff) to extract features utilizing the same pipeline as \code{\link{axDistmap}} or directly utilizes precomputed feature data (.txt).
#'   \item \strong{SVM classification}: Utilizes the loaded SVM model to classify axonal regions into "Degenerate" or "Intact" based on morphological features.
#'   \item \strong{Index calculation}: Calculates Axon Integrity Index (AII) and Degeneration Index (DI) based on classified areas:
#'         \deqn{AII = \frac{\text{Area of Intact Axons}}{\text{Total Axonal Area}}}
#'         \deqn{DI = \frac{\text{Area of Degenerate Axons}}{\text{Total Axonal Area}}}
#'   \item \strong{Result saving}: Exports summary and optional detailed prediction files with unique filenames to prevent overwriting.
#' }
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
#'   \item The function requires a pre-trained SVM model compatible with the feature set used by \code{\link{axDistmap}} and \code{\link{axSvm}}.
#'   \item Required features for SVM classification: m.eccentricity, s.radius.sd, h.sva.s2, h.idm.s1, h.sen.s1, m.majoraxis.
#'   \item When \code{imprtImg = TRUE}, images are processed using the same pipeline as \code{\link{axDistmap}} with the specified \code{resizeW} and \code{subBack} parameters.
#'   \item Single image prediction results are exported only if \code{expSip = TRUE}.
#'   \item Compatible with feature data and models from \code{\link{axDistmap}} and \code{\link{axSvm}}.
#'   \item Output files include timestamp in ISO 8601 format (YYYY-MM-DD).
#' }
#'
#' @examples
#' \dontrun{
#' # Interactive mode: process .tiff images with GUI dialogs
#' # NOTE: This example requires a GUI environment for interactive folder selection.
#' axQnt(imprtImg = TRUE, expSip = TRUE)
#'
#'
#' # Utilize package-included image folder and output to temporary directory
#' img_dir <- system.file("extdata", "Degenerate_Images", package = "AiES")
#' svm_model <- system.file("extdata", "svm_example_model.svm", package = "AiES")
#' axQnt(imprtImg = TRUE, svm_model_path = svm_model, input_dirs = img_dir, output_dir = tempdir())
#'
#' # Process with custom image resize and background threshold
#' img_dir <- system.file("extdata", "Intact_Images", package = "AiES")
#' svm_model <- system.file("extdata", "svm_example_model.svm", package = "AiES")
#' axQnt(imprtImg = TRUE, resizeW = 700, subBack = 50,
#' svm_model_path = svm_model, input_dirs = img_dir, output_dir = tempdir())
#' }
#'
#' # Utilize package-included precomputed feature data from .txt files
#' # and output to temporary directory
#' txt_dir <- system.file("extdata", "Degenerate_txt", package = "AiES")
#' svm_model <- system.file("extdata", "svm_example_model.svm", package = "AiES")
#' axQnt(imprtImg = FALSE, svm_model_path = svm_model, input_dirs = txt_dir, output_dir = tempdir())
#'
#'
#' @import stringr dplyr
#' @importFrom data.table fread fwrite rbindlist data.table
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment
#' @importFrom stats predict na.omit
#' @export
#'
#'

######################Codes for test run
##library("data.table")
##library("dplyr")
##library("stringr")
##library("EBImage")
##library("stats")

##以下は不要250523
##library("ggpubr")
##library("colorspace")

##subBack <- 30
##imprtImg = FALSE
##expSip = TRUE
###############Definition of variable######################


axQnt <- function(imprtImg = TRUE,
                  subBack = 30,
                  resizeW = 900,
                  binaryImg = FALSE,
                  expSip = TRUE,
                  svm_model_path = NULL,
                  input_dirs = NULL,
                  output_dir = NULL){

    sdate <- Sys.Date()
    txtGroup <- c("Degenerate","Intact")
    #Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
    ftrSvm <- c("m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
    rltName <- c("FileName","AxonIntegrityIndex","DegenerationIndex")

    ###############Definition of function######################
    ##Function 1
    # クロスプラットフォームなファイル・フォルダ選択
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

    ##Function 2
    select_files <- function(caption = "Select files", pattern = NULL) {
        if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.files(caption = caption, multi = TRUE))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_open(title = caption, multiple = TRUE)$res)
        } else {
            cat("Enter file path(s), separated by comma: ")
            return(strsplit(readline(), ",")[[1]])
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

    # 1. SVMモデルのロード
    if (is.null(svm_model_path)) {
        svm_model_path <- select_files("Select a SVM model file")
        if (length(svm_model_path) == 0 || is.na(svm_model_path)) return(message("SVM model not selected."))
    }

    loaded_objects <- load(svm_model_path)
    svm_model <- NULL##変数初期化
    # ロードされたオブジェクトからSVMモデルを探す
    for (obj_name in loaded_objects) {
        obj <- get(obj_name)
        if (inherits(obj, "svm")) {
            svm_model <- obj
            break
        }
    }

    if (is.null(svm_model)) {
        stop("No SVM model found in the loaded file", svm_model_path)
    }

    # 2. 入力ディレクトリの選択
    if (is.null(input_dirs)) {
        input_dirs <- character(0) # 初期化
        continue_selection <- "y"
        repeat {
            selected_folder <- select_folder(if (imprtImg) "Select folder containing TIFF images" else "Select folder containing TXT feature files")
            if (is.null(selected_folder)) {
                message("Folder selection cancelled.")
                break # ループを抜ける
            }
            # 選択されたフォルダが既にリストにないか確認
            if (!(selected_folder %in% input_dirs)) {
                input_dirs <- c(input_dirs, selected_folder)
                message(sprintf("Added folder: %s", selected_folder))
            } else {
                message(sprintf("Folder already selected: %s. Skipping.", selected_folder))
            }
            # フォルダ選択を続けるか確認
            continue_selection <- readline(prompt = "Select another folder? (y/n): ")
            if (tolower(continue_selection) != "y") {  #大文字を小文字に変換
                break # ループを抜ける
            }
        }
        # フォルダが一つも選択されなかった場合
        if (length(input_dirs) == 0) {
            message("No input folder selected. Exiting function.")
            return(NULL)
        }
    }

    # 3. 出力ディレクトリの選択
    if (is.null(output_dir)) {
        output_dir <- select_folder("Select output folder")
        if (is.null(output_dir) || output_dir == "") output_dir <- getwd()
    }
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # 4. ファイルリスト取得（複数フォルダ＆サブフォルダ対応）
    if (imprtImg) {
        input_files <- unlist(lapply(input_dirs, function(dir) {
            list.files(dir, pattern = "\\.tiff?$", full.names = TRUE, recursive = TRUE)
        }))
    } else {
        input_files <- unlist(lapply(input_dirs, function(dir) {
            list.files(dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
        }))
    }
   if (length(input_files) == 0) return(message("No input files found."))


    # 5. 画像または特徴量テーブルの処理
    results_summary <- data.table::data.table(FileName = character(0), AxonIntegrityIndex = numeric(0), DegenerationIndex = numeric(0))

    # lapplyで画像または特徴量テーブルの処理
    results_list <- lapply(input_files, function(file) {
        # 画像処理または特徴量読み込み
        if (imprtImg) {
            message(sprintf("Processing file: %s", file))
            tmpImage <- EBImage::readImage(file)
            tmpImage <- EBImage::resize(tmpImage, w = resizeW)
            contrastST <- (1/(1+(0.5/tmpImage[,,1])^5))
            medFltr <- EBImage::medianFilter(contrastST,1)
            logTrsf <- (log1p(medFltr)/log1p(max(medFltr)))
            logTrsfLog <- (log1p(logTrsf)/log1p(max(logTrsf)))
            thrWBInv <- EBImage::thresh(logTrsfLog, 8, 8, -0.04)
            kern <- EBImage::makeBrush(1, shape="box")
            mrphOprt <- 1 - EBImage::opening(EBImage::closing(thrWBInv, kern), kern)
            dm <- EBImage::distmap(mrphOprt)
            nDm <- EBImage::normalize(dm)
            bnrySeg <- EBImage::bwlabel(nDm)
            sdat <- as.data.frame(EBImage::computeFeatures.shape(bnrySeg))
            rma <- subBack
            rmNum <- which(sdat$s.area <= rma)
            bnrySeg <- EBImage::rmObjects(bnrySeg, rmNum)
            sdat <- as.data.frame(EBImage::computeFeatures.shape(bnrySeg))
            hdat <- as.data.frame(EBImage::computeFeatures.haralick(bnrySeg,nDm))
            mdat <- as.data.frame(EBImage::computeFeatures.moment(bnrySeg))
            singleData <- data.frame(cbind(sdat, hdat, mdat))
            singleData <- na.omit(singleData)
            singleData <- singleData[, c("s.area", ftrSvm), drop = FALSE]
            # 必要に応じてファイル保存
        } else {
            singleData <- tryCatch(data.table::fread(file), error = function(e) NULL)
            if (is.null(singleData)) return(NULL)
        }
        # SVM分類
        if (!all(ftrSvm %in% colnames(singleData))) {
            message("Required features missing in file: ", file)
            return(NULL)
        }
        ## error 対策
        singleData <- as.data.frame(singleData)
        #pred <- predict(svm_model, singleData[, ftrSvm, with = FALSE], type = "class") #data.table
        #pred <- predict(svm_model, singleData[, ..ftrSvm], type = "class") #data.frame
        pred <- predict(svm_model, singleData[, ftrSvm, drop = FALSE], type = "class") #data.frame

        # AII/DI計算
        intact_area <- sum((pred == "Intact") * singleData$s.area)
        deg_area <- sum((pred == "Degenerate") * singleData$s.area)
        total_area <- intact_area + deg_area
        AII <- if (total_area > 0) intact_area / total_area else NA
        DI <- if (total_area > 0) deg_area / total_area else NA
        results_summary <- rbind(results_summary, data.table::data.table(FileName = basename(file), AxonIntegrityIndex = AII, DegenerationIndex = DI))
        # 単画像予測結果の保存
        if (expSip) {
            sip_outfile <- file.path(output_dir, paste0(basename(file), "_", sdate, "_SIP.csv"))
            sip_outfile <- generate_unique_filename(sip_outfile)
            data.table::fwrite(data.frame(Pred = pred, singleData), file = sip_outfile, sep = ",")
        }
        return(data.table(FileName = basename(file), AxonIntegrityIndex = AII, DegenerationIndex = DI))
    })

    # リストをdata.tableにまとめる
    results_summary <- data.table::rbindlist(results_list, fill = TRUE)

    # サマリー出力
    summary_outfile <- file.path(output_dir, paste0("Summary_", sdate, ".csv"))
    summary_outfile <- generate_unique_filename(summary_outfile)
    data.table::fwrite(results_summary, file = summary_outfile, sep = ",")
    message("Summary file saved: ", summary_outfile)
    invisible(results_summary)


}#####end of function
