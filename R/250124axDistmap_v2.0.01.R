# Copyright 2025 Your Company Name
# BSD 3-Clause License (see LICENSE file)
#' @title Create distance map and binary image from TIFF image file
#' @description
#' `axDistmap` processes TIFF image files to create distance maps and optionally binary images.
#' It also computes various image features and exports them as a text file.
#'
#' @param subBack Numeric. Size of background objects to subtract (default: 30 pixels).
#' @param Binary Logical. If TRUE, exports binary image files (default: FALSE).
#' @param allFeatures Logical. If TRUE, exports data of all computed features (default: FALSE).
#' @param imgType Character. Output image format: "png", "jpg", or "tiff" (default: "tiff").
#' @param folder_paths Character vector. A character vector of folder paths to process. If NULL, a folder selection dialog will be shown. (default: "NULL").
#'
#' @return This function doesn't return a value directly, but produces the following outputs:
#' \itemize{
#'   \item A distance map image file (format specified by `imgType`)
#'   \item A text file containing computed image features (named "'original_filename'_'current_date'_ImageData.txt")
#'   \item (Optional) A binary image file if `Binary = TRUE` (format specified by `imgType`)
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Reads and resizes the input TIFF image
#' 2. Applies various image processing techniques (contrast adjustment, filtering, thresholding)
#' 3. Computes a distance map
#' 4. Computes shape, Haralick, and moment features
#' 5. Exports the processed images as image files
#' 6. Exports the computed features as a tab-separated text file
#'
#' @section Feature Export:
#' The function always exports computed features as a text file.
#' If allFeatures = FALSE, it exports a subset of features including:
#' - s.area (EBImage::computeFeatures.shape)
#' - m.eccentricity (EBImage::computeFeatures.moment)
#' - s.radius.sd (EBImage::computeFeatures.shape)
#' - h.sva.s2 (EBImage::computeFeatures.haralick)
#' - h.idm.s1 (EBImage::computeFeatures.haralick)
#' - h.sen.s1 (EBImage::computeFeatures.haralick)
#' - m.majoraxis (EBImage::computeFeatures.moment)
#'
#' If allFeatures = TRUE, it exports all computed features from EBImage's computeFeatures
#' functions (shape, moment, and haralick) plus an additional 'Cir' feature.
#' The 'Cir' feature represents Circularity and is calculated as:
#' Cir = (s.area * pi * 4) / (s.perimeter^2), where s.area and s.perimeter are from
#' EBImage::computeFeatures.shape.
#'
#' @note
#' - The function will prompt the user to select a directory containing TIFF files.
#' - It processes all TIFF files in the selected directory.
#' - Once there are no unprocessed files left in the selected directory, the function will prompt
#'   the user to choose whether to process another folder.
#' - Processing will continue until the user cancels the operation.
#' - Output files (images and feature data) are saved in the same directory as the input files.
#' - The feature data text file is named using the format: "'original_filename'_'current_date'_ImageData.txt"
#' - If allFeatures = FALSE, only features required for the support vector machine learning in this package are exported.
#' - If allFeatures = TRUE, all features from EBImage package plus Circularity are exported.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' # Exports only features needed for SVM learning
#' axDistmap()
#'
#' # Create binary images and export all EBImage features plus Circularity as PNG
#' axDistmap(subBack = 50, Binary = TRUE, allFeatures = TRUE, imgType = "png")
#'
#' # Process images and export as TIFF without binary images
#' # Only exports features needed for SVM learning
#' axDistmap(subBack = 20, Binary = FALSE, allFeatures = FALSE, imgType = "tiff")
#'
#' # Export all EBImage features plus Circularity without creating binary images
#' axDistmap(Binary = FALSE, allFeatures = TRUE)
#' }
#'
#' @import stringr colorspace dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage
#' @importFrom utils write.table choose.files
#' @importFrom ggpubr mutate
#' @importFrom stats na.omit
#' @importFrom fs dir_ls
#'
#' @export


axDistmap <- function(subBack = 30, Binary = FALSE, allFeatures = FALSE, imgType = "tiff", folder_paths = NULL){

    ###imgType check
    if(!(imgType %in% c("tiff","png","jpg"))) imgType <- "tiff"

    #### #test code for debug
    ##library("EBImage")
    ##library("stringr")
    ##library("colorspace")
    ##library("EBImage")
    ##library("dplyr")
    ##subBack <- 30
    ##Binary = "TRUE"
    ##allFeatures = "TRUE"
    ##imgType = "tiff"
    ##folder_paths = c("D:/NCNP/01 R&D Eng/01 Projects/04 R package/02 R&D/RStudio/AiES/inst,
    ##                 "D:/NCNP/01 R&D Eng/01 Projects/04 R package/02 R&D/RStudio/AiES/inst/")
    #########################
    ###############Definition of variable######################

    sdate <- Sys.Date()
    ###############Definition of function######################
    ##Function 1
    ## Functions are executed depending on the user environment
    select_folder <- function() {
        if (rstudioapi::isAvailable()) {
            return(rstudioapi::selectDirectory(caption = "Select folder"))
        } else if (requireNamespace("tcltk", quietly = TRUE)) {
            return(tcltk::tk_choose.dir(caption = "Select folder"))
        } else if (requireNamespace("svDialogs", quietly = TRUE)) {
            return(svDialogs::dlg_dir(title = "Select folder")$res)
        } else {
            message("No GUI folder selection method available")
            return(NULL)
        }
    }

    # フォルダの選択または指定
    if (is.null(folder_paths)) {
        folder_paths <- select_folder()
        if (is.null(folder_paths)) {
            message("No folder selected. Exiting function.")
            return(NULL)
        }
    } else {
        folder_paths <- as.list(folder_paths)
    }


    # 各フォルダを処理
    results <- lapply(folder_paths, function(folder_path) {
        if (!dir.exists(folder_path)) {
            message(sprintf("Specified folder does not exist: %s. Skipping.", folder_path))
            return(NULL)
        }

        # TIFFファイル(.tiff or .tif)の取得
        tiff_files <- fs::dir_ls(folder_path, regexp = "\\.tiff?$", recurse = TRUE)

        if (length(tiff_files) == 0) {
            message(sprintf("No TIFF files found in the folder: %s. Skipping.", folder_path))
            return(NULL)
        }

        # 複数のファイルを処理
        process_file <- function(file) {
            message(sprintf("Processing file: %s", file))


            tmpImage <- readImage(file)#file_list.name

                # 二値化画像かどうかをチェック
            if (all(tmpImage %in% c(0, 1))) {
              message(sprintf("Skipping binary image: %s", file))
              return(NULL)
            }

                # グレースケール変換
            if (length(dim(tmpImage)) == 3) {
              tmpImage <- channel(tmpImage, "gray")
            }


            tmpImage <- resize(tmpImage, w = 900)#15Jun11:696 22Apr15:900
            contrastST <- (1/(1+(0.5/tmpImage[,,1])^5))
            medFltr <- medianFilter(contrastST,1)
            logTrsf <- (log1p(medFltr)/log1p(max(medFltr)))
            logTrsfLog <- (log1p(logTrsf)/log1p(max(logTrsf)))
            thrWBInv <- 1 - thresh(logTrsfLog, 8, 8, -0.04)
            kern <- makeBrush(1, shape="box")
            mrphOprt <- opening(closing(thrWBInv, kern), kern)
            dm <- distmap(mrphOprt)
            nDm <- normalize(dm)
            bnrySeg <- bwlabel(nDm)
            sdat <- as.data.frame(computeFeatures.shape(bnrySeg))
            #######Threshold for eliminating small objects that cannot be classified by image processing
            rma <- subBack
            rmNum <- which(sdat$s.area <= rma)
            bnrySeg <- rmObjects(bnrySeg, rmNum)
            sdat <- as.data.frame(computeFeatures.shape(bnrySeg)) #option with
            hdat <- as.data.frame(computeFeatures.haralick(bnrySeg,nDm))
            mdat <- as.data.frame(computeFeatures.moment(bnrySeg))
            singleData <- data.frame(cbind(sdat, hdat, mdat))
            invisible({rm(list=c("sdat", "hdat", "mdat"));gc();gc()})
            singleData <- na.omit(singleData)
            #if (is.numeric(singleData[,1]) == FALSE){
            #  singleData <- sapply(singleData,num_f)
            #}

            if(allFeatures == FALSE){
                singleData <- data.frame(s.area = singleData$s.area,
                                         m.eccentricity = singleData$m.eccentricity,
                                         s.radius.sd = singleData$s.radius.sd,
                                         h.sva.s2 = singleData$h.sva.s2,
                                         h.idm.s1 = singleData$h.idm.s1,
                                         h.sen.s1 = singleData$h.sen.s1,
                                         m.majoraxis = singleData$m.majoraxis,
                                         stringsAsFactors = TRUE)
            }else {
                singleData <- data.frame(cbind(singleData,Cir = singleData$s.area*pi*4/singleData$s.perimeter^2) , stringsAsFactors = TRUE)
            }

            #####export data file
            write.table(singleData, paste0(file,"_",sdate,"_ImageData.txt"),
                        sep="\t",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
            message(sprintf("%s_%s_ImageData.txt",file,sdate))

            if (Binary == TRUE){
                dmrmabw <- 1*(bnrySeg > 0)# %>% # binary image
                file.name <- paste0(file,"_",sdate,  "_Binary.", imgType)
                writeImage(dmrmabw, file.name,type = imgType,quality = 100)
                message(file.name)

            }
            nDmrmadm <- nDm*(bnrySeg > 0)# %>% # distance map

            file.name <- paste0(file,"_",sdate,  "_DistMap.", imgType)
            writeImage(nDmrmadm, file.name,type = imgType,quality = 100)
            message(file.name)
        }


        folder_results <- lapply(tiff_files, process_file)

        message(sprintf("All files in folder %s processed successfully.", folder_path))

        return(folder_results)
    })

    # 処理完了メッセージ
    message("All folders processed successfully.")

    ##return(NULL)
}














