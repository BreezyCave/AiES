# Copyright 2025 National Center of Neurology and Psychiatry
# BSD 3-Clause License (see LICENSE file)
#' @title Axon integrity index quantification
#' @description
#' `axQnt` quantifies the axonal integrity index (Axon Integrity Index, AII) and degeneration index (Degeneration Index, DI)
#' from either image files (.tiff) or precomputed data files (.txt). The function uses a pre-trained SVM model to classify
#' axonal regions and calculates the indices based on the classified areas. It also exports processed data and summary files.
#'
#'
#' @import stringr dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment
#' @importFrom stats predict na.omit
#' @importFrom ggpubr mutate
#' @importFrom utils write.table choose.files read.table
#'
#' @param imprtImg Logical. If TRUE, imports image files (.tiff) for processing. If FALSE, imports precomputed data files (.txt) (default: TRUE).
#' @param subBack Numeric. Threshold for subtracting background objects based on size (default: 30 pixels).
#' @param expSip Logical. If TRUE, exports single image prediction results as .csv files (default: TRUE).
#'
#' @return The function does not return values directly but generates the following outputs:
#' \itemize{
#'   \item A summary .csv file containing the calculated Axon Integrity Index (AII) and Degeneration Index (DI) for each file.
#'   \item (Optional) Single image prediction results as .csv files if `expSip = TRUE`.
#'   \item (Optional) Processed image data as .txt files if `imprtImg = TRUE`.
#' }
#'
#' @details
#' The function processes input data as follows:
#' 1. Prompts the user to select a pre-trained SVM model file.
#' 2. Depending on `imprtImg`, it either processes image files (.tiff) to extract features or directly uses precomputed feature data (.txt).
#' 3. Uses the loaded SVM model to classify axonal regions into "Degenerate" or "Intact".
#' 4. Calculates Axon Integrity Index (AII) and Degeneration Index (DI) based on classified areas:
#'    \[
#' \deqn{AII = \frac{\text{Area of Intact Axons}}{\text{Total Axonal Area}}}
#' \deqn{DI = \frac{\text{Area of Degenerate Axons}}{\text{Total Axonal Area}}}
#'    \]
#' 5. Exports results as summary and optional detailed prediction files.
#'
#' @note
#' - The function uses a recursive process to allow continuous processing of multiple directories until the user cancels.
#' - Processed image data is exported only if `imprtImg = TRUE`.
#' - Single image prediction results are exported only if `expSip = TRUE`.
#'
#' @examples
#' \dontrun{
#' # Process .tiff images to calculate Axon Integrity Index and export single image predictions
#' axQnt(imprtImg = TRUE, subBack = 30, expSip = TRUE)
#'
#' # Use precomputed feature data from .txt files without exporting single image predictions
#' axQnt(imprtImg = FALSE, subBack = 20, expSip = FALSE)
#'
#' # Adjust background subtraction threshold and process images
#' axQnt(imprtImg = TRUE, subBack = 50)
#' }
#'
#' @export
#'


axQnt <- function(imprtImg = TRUE, subBack = 30, expSip = TRUE){

    ######################Codes for test run
    ##library("dplyr")
    ##library("stats")
    ##library("stringr")
    ##library("ggpubr")
    ##library("EBImage")
    ##library("stringr")
    ##library("colorspace")
    ##library("EBImage")
    ##library("dplyr")
    ##subBack <- 30
    ##imprtImg = FALSE
    ##expSip = TRUE
    ###############Definition of variable######################

    txtGroup <- c("Degenerate","Intact")
    #Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
    ftrSvm <- c("m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
    rltName <- c("FileName","AxonIntegrityIndex","DegenerationIndex")
    sdate <- Sys.Date()

    ###############Definition of function######################
    ##Function 1
    createEmptyDF <- function( nrow, ncol, colnames = c() ){
        if( missing( ncol ) && length( colnames ) > 0 ){
            ncol <- length( colnames )
        }
        data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
    }

    ##Function 2
    isTxt <- function(x) regexpr('\\.txt$', x) > 0
    ##Function 3
    isTif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0
    ##Function 4
    #num_f <- function(x){
    #  x <- as.numeric(levels(x))[x]
    #}

    ##Function 5
    ####function to extract dirInfo
    dirInfo <- function(x){
        str_replace_all(x, pattern = "\\\\", replacement="/") %>%
            str_split( "/") %>%
            unlist() %>%
            as.data.frame()  %>%
            mutate(level = row_number())
    }

    ##Function 6
    ########Processing Image
    procImage <- function(x){
        tmpImage <- readImage(x)#file_list.name
        tmpImage <- resize(tmpImage, w = 900)#15Jun11:696 22Apr15:900
        contrastST <- (1/(1+(0.5/tmpImage[,,1])^5))
        medFltr <- medianFilter(contrastST,1)
        logTrsf <- (log1p(medFltr)/log1p(max(medFltr)))
        logTrsfLog <- (log1p(logTrsf)/log1p(max(logTrsf)))
        thrWBInv <- thresh(logTrsfLog, 8, 8, -0.04)
        kern <- makeBrush(1, shape="box")
        mrphOprt <- 1 - opening(closing(thrWBInv, kern), kern)
        dm <- distmap(mrphOprt)
        nDm <- normalize(dm)
        bnrySeg <- bwlabel(nDm)
        sdat <- as.data.frame(computeFeatures.shape(bnrySeg))
        #######Threshold for eliminating small objects that cannot be classified by image processing
        rma <- subBack
        rmNum <- which(sdat$s.area <= rma)
        bnrySeg <- rmObjects(bnrySeg, rmNum)
        sdat <- as.data.frame(computeFeatures.shape(bnrySeg))
        hdat <- as.data.frame(computeFeatures.haralick(bnrySeg,nDm))
        mdat <- as.data.frame(computeFeatures.moment(bnrySeg))

        singleData <- data.frame(cbind(sdat, hdat, mdat))
        invisible({rm(list=c("sdat", "hdat", "mdat"));gc();gc()})
        singleData <- na.omit(singleData)
        #if (is.numeric(singleData[,1]) == FALSE){
        #  singleData <- sapply(singleData,num_f)
        #}

        singleData <- data.frame(s.area = singleData$s.area,
            m.eccentricity = singleData$m.eccentricity,
            s.radius.sd = singleData$s.radius.sd,
            h.sva.s2 = singleData$h.sva.s2,
            h.idm.s1 = singleData$h.idm.s1,
            h.sen.s1 = singleData$h.sen.s1,
            m.majoraxis = singleData$m.majoraxis,
            stringsAsFactors = TRUE)

        #####export data file
        write.table(singleData, paste0(x,sdate,"_ImageData.txt"),
                    sep="\t",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
        message(sprintf("%s_%s_ImageData.txt",x,sdate))
    }

    ##Function 7
    ########Processing Text data
    procTxt <- function(x, y){
        Data_shq <- read.table(x, header=TRUE, sep="\t")
        if(! FALSE %in% (ftrSvm %in% colnames(Data_shq))){
            Pred <- svmModelLoaded %>%
                get() %>%
                predict(Data_shq, type="class", probability = FALSE) %>%
                as.character()
            if(expSip == TRUE){
                write.table(data.frame(Pred = Pred,Data_shq),
                            paste0(x, "_", sdate,"_SIP.csv"),
                            sep=",",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
                message(sprintf( "%s_%s_SIP.csv", x,sdate))
            }
            AxII <- sum((Pred %in% txtGroup[2])*Data_shq$s.area)/(sum((Pred %in% txtGroup[2])*Data_shq$s.area)+sum((Pred %in% txtGroup[1])*Data_shq$s.area))
            DegI <- sum((Pred %in% txtGroup[1])*Data_shq$s.area)/(sum((Pred %in% txtGroup[2])*Data_shq$s.area)+sum((Pred %in% txtGroup[1])*Data_shq$s.area))
            y[nrow(y)+1,] <- c(x, AxII, DegI)
            return(y)
        }else{
            message(x, ": This file does not contain the required data \n")
        }
    }

    ##Function 8
    ######Selecting the Directory#######
    procFiles <- function() {  ##recursive function
        # file full path
        ifelse(imprtImg == TRUE,
            filePath2 <- choose.files(caption = "Select any tiff file to set the directory",
                                         multi=FALSE),
            filePath2 <- choose.files(caption = "Select any text file to set the directory",
                                         multi=FALSE))
        if(length(filePath2)==0) return(("File not selected.\n"))

        # extract directory info from the selected file
        dirInfo2 <- dirInfo(filePath2)
        colnames(dirInfo2) <- c("Dir_Name", "Level")
        setwd(str_c(dirInfo2$Dir_Name[seq_len(nrow(dirInfo2)-1)], collapse = "/"))


        if(imprtImg == TRUE){
            ################Calculate axon integrity index from .tiff data#####
            lapply(list.files()[isTif(list.files())], procImage)
        }

        ##
        ################Calculate axon integrity index from .txt data#####
        rltSummary <- createEmptyDF(0, colnames = rltName)
        rltSummary <- lapply(list.files()[isTxt(list.files())], procTxt, y=rltSummary)
        rltSummary <- do.call(rbind, rltSummary)

        ##################Export the Summary data (csv)##############################
        write.table(rltSummary, paste0(dirInfo2$Dir_Name[(nrow(dirInfo2)-1)], "_", sdate,"_Summary.csv"),
                    sep=",",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
        message(sprintf("%s_%s_Summary.csv", dirInfo2$Dir_Name[(nrow(dirInfo2)-1)],sdate))


        procFiles()##recursive function
    }

    ################################
    ##Excecute
    #############import svm model #######
    ######Selecting the Directory1#######
    # file full path
    filePath1 <- choose.files(caption = "Select a SVM model",
        multi=FALSE)
    if (length(filePath1)==0) {
        return(message("File not selected.\n"))##Error countermeasure code
    }else{
        tryCatch({
            dirInfo1 <- dirInfo(filePath1)
            colnames(dirInfo1) <- c("Dir_Name", "Level")
            svmModelLoaded <- load(filePath1)
        }, error = function(e){
            stop(" Error in load file '",
                dirInfo1$Dir_Name[nrow(dirInfo1)],"'\n This file dose not have SVM model. \n") ##* NOTE:  Avoid redundant 'stop' and 'warn*' in signal conditions
        })
    }

    #####################################
    procFiles()



}#####end of function








