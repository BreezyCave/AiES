# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#' @title Build SVM Classifier for Axon Image Feature Classification
#' @description
#' Trains a Support Vector Machine (SVM) classifier to distinguish between degenerative
#' and intact axon states using image feature data. Supports both direct path input
#' and interactive GUI selection.
#'
#' @param degenerate_paths Character vector. Path(s) to degenerative feature data
#'        (folder containing .txt files or direct file paths). If NULL, GUI selection
#'        dialog will be shown. Default: NULL
#' @param intact_paths Character vector. Path(s) to intact feature data
#'        (folder containing .txt files or direct file paths). If NULL, GUI selection
#'        dialog will be shown. Default: NULL
#' @param output_data_path Character. Path to save combined feature data.
#'        If NULL, GUI save dialog will be shown. Default: NULL
#' @param output_model_path Character. Path to save trained SVM model.
#'        If NULL, GUI save dialog will be shown. Default: NULL
#' @param nCst Numeric. SVM cost parameter (C-value) controlling margin hardness.
#'        Higher values increase model complexity. Default: 3
#' @param nGmm Numeric. SVM gamma parameter (γ-value) controlling RBF kernel width.
#'        Smaller values mean larger kernel radius. Default: 0.1
#' @param nCrss Integer. Number of folds for cross-validation.
#'        Recommended values 5-10. Default: 5
#'
#' @return This function doesn't return values directly, but produces the following outputs:
#' \itemize{
#'   \item Combined feature data file for machine learning (.txt format).
#'   \item Trained SVM model file (.svm format)
#' }
#'
#' @details
#' The function implements a complete machine learning pipeline:
#' \enumerate{
#'   \item \strong{Input Handling}: Accepts both directory paths and direct file paths
#'   \item \strong{Data Loading}: Uses data.table::fread() for efficient large file handling
#'   \item \strong{Feature Selection}: Focuses on 7 key morphological features:
#'         \itemize{
#'           \item m.eccentricity (elliptical eccentricity)
#'           \item s.radius.sd (radial distribution uniformity)
#'           \item h.sva.s2 (Sum Variance scale=2)
#'           \item h.idm.s1 (local homogeneity scale=1)
#'           \item h.sen.s1 (structural complexity scale=1)
#'           \item m.majoraxis (length)
#'           \item s.area (area size)
#'         }
#'   \item \strong{Model Training}: Utilizes radial basis function (RBF) kernel SVM
#'   \item \strong{Model Evaluation}: Includes built-in k-fold cross-validation
#' }
#'
#' @section Feature Selection:
#' The function uses the following features for SVM training:
#' "Group", "m.eccentricity", "s.radius.sd", "h.sva.s2", "h.idm.s1", "h.sen.s1", "m.majoraxis"
#'
#' @note
#' \itemize{
#'   \item Minimum 8GB RAM recommended for typical datasets
#'   \item Compatible with feature data from \code{\link{axDistmap}}
#'   \item Output files include timestamp in ISO 8601 format (YYYY-MM-DD)
#'   \item Model files can be reloaded using base::load()
#' }
#'
#' @examples
#' \dontrun{
#' # Interactive mode with GUI prompts
#' axSvm()
#'
#' # Direct path specification
#' axSvm(
#'   degenerate_paths = "data/deg_features",
#'   intact_paths = "data/intact_features",
#'   output_data_path = "results/features_2025.tsv",
#'   output_model_path = "models/svm_model_2025.svm"
#' )
#'
#' # Custom hyperparameters
#' axSvm(nCst = 5, nGmm = 0.05, nCrss = 10)
#' }
#'
#' @importFrom data.table fread rbindlist
#' @importFrom e1071 svm
#' @importFrom utils choose.files
#' @export




axSvm <- function(nCst = 3, nGmm = 0.1, nCrss=5,
                  degenerate_paths = NULL,
                  intact_paths = NULL,
                  output_data_path = NULL,
                  output_model_path = NULL
){

    ######################Codes for test run
    ##library("dplyr")
    ##library("e1071")
    ##library("stringr")
    ##library("ggpubr")
    ##Image = TRUE
    ##nCst = 3
    ##nGmm = 0.1
    ##nCrss=5
    #######################################


    ##############Definition of variable######################

    sdate <- Sys.Date()
    ftrTmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
    txtGroup <- c("Degenerate","Intact")
    #mainFldr <- if(exists("mainFldr")) mainFldr else c(NA_character_, NA_character_)

    ###############Definition of function######################
    ##Function 1 フォルダ選択
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

    ##Function 2 ファイル選択
    select_files <- function(caption = "Select files") {
        if (.Platform$OS.type == "windows") {
            return(utils::choose.files(caption = caption, multi = TRUE))
        } else if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.files(caption = caption, multi = TRUE))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_open(title = caption, multiple = TRUE)$res)
        } else {
            stop("No GUI file selection method available on this OS.")
        }
    }

    # 1. 入力データの取得
    if (is.null(degenerate_paths)) {
        message("Select folder or .txt files for Degenerate condition")
        degenerate_paths <- select_folder("Select Degenerate folder")
        if (is.na(degenerate_paths) || is.null(degenerate_paths)) return(message("No folder selected."))
        degenerate_files <- list.files(degenerate_paths, pattern = "\\.txt$", full.names = TRUE)
    } else if (dir.exists(degenerate_paths)) {
        degenerate_files <- list.files(degenerate_paths, pattern = "\\.txt$", full.names = TRUE)
    } else {
        degenerate_files <- degenerate_paths
    }
    if (length(degenerate_files) == 0) return(message("No degenerate files found."))

    if (is.null(intact_paths)) {
        message("Select folder or .txt files for Intact condition")
        intact_paths <- select_folder("Select Intact folder")
        if (is.na(intact_paths) || is.null(intact_paths)) return(message("No folder selected."))
        intact_files <- list.files(intact_paths, pattern = "\\.txt$", full.names = TRUE)
    } else if (dir.exists(intact_paths)) {
        intact_files <- list.files(intact_paths, pattern = "\\.txt$", full.names = TRUE)
    } else {
        intact_files <- intact_paths
    }
    if (length(intact_files) == 0) return(message("No intact files found."))


    # 2. データ読み込み（大規模対応）
    read_data <- function(files, group) {
        datalist <- lapply(files, function(f) {
            # freadはdata.tableの高速読み込み
            dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
            if (is.null(dt)) return(NULL)
            dt$Group <- group
            dt
        })
        datalist <- datalist[!sapply(datalist, is.null)]
        if (length(datalist) == 0) return(NULL)
        data.table::rbindlist(datalist, fill = TRUE)
    }
    deg_data <- read_data(degenerate_files, "Degenerate")
    intact_data <- read_data(intact_files, "Intact")
    if (is.null(deg_data) || is.null(intact_data)) return(message("Failed to read data."))



    # 3. データ結合・特徴量選択
    dataSvm <- rbind(deg_data, intact_data)
    dataSvm <- dataSvm[, ..ftrTmp]
    dataSvm$Group <- as.factor(dataSvm$Group)
    rm(deg_data, intact_data); gc()

    # 4. SVMモデル構築
    svmModel <- e1071::svm(
        Group ~ .,
        data = dataSvm,
        method = "C-classification",
        probability = TRUE,
        kernel = "radial",
        gamma = nGmm,
        cost = nCst,
        cross = nCrss
    )

    # 5. 出力ファイルの指定
    if (is.null(output_data_path)) {
        default_name <- paste0(sdate, "_Extracted_data_for_ML.txt")
        output_data_path <- utils::choose.files(default = default_name, caption = "Save extracted data file", multi = FALSE)
        if (length(output_data_path) == 0 || output_data_path == "") return(message("Canceled"))
    }
    data.table::fwrite(dataSvm, file = output_data_path, sep = "\t")
    message("Extracted data saved: ", output_data_path)

    if (is.null(output_model_path)) {
        default_name <- paste0(sdate, "_AxClassifier.svm")
        output_model_path <- utils::choose.files(default = default_name, caption = "Save SVM model", multi = FALSE)
        if (length(output_model_path) == 0 || output_model_path == "") return(message("Canceled"))
    }
    save(svmModel, file = output_model_path)
    message("SVM model saved: ", output_model_path)

    invisible(svmModel)
}###end of function









