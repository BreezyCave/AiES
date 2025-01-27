# Copyright 2024 Your Company Name
# BSD 3-Clause License (see LICENSE file)
#' @title creating distance map and binary image from tiff image file
#' @description \code{axDistmap} create images of distance map and binary image (option) from tiff image files
#' @import stringr
#' @import colorspace
#' @import dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage
#' @importFrom utils write.table
#' @importFrom ggpubr mutate
#' @importFrom utils choose.files
#' @importFrom stats na.omit
#' @param subBack 30: subtract background objects (default 30 pixels)
#' @param Binary TRUE: exporting binary image files
#' @param allFeatures TRUE: exporting data of all features
#' @param imgType png, jpg, tiff
#' @return return the image of distancemap and data of features
#' @export
#' @examples
#' # axDistmap(subBack = 30, Binary = FALSE, allFeatures = FALSE, imgType = png)


axDistmap <- function(subBack = 30, Binary = FALSE, allFeatures = FALSE, imgType = "png"){

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
    #########################
    ###############Definition of variable######################

    sdate <- Sys.Date()
    ###############Definition of function######################
    ##Function 1
    #num_f <- function(x){
    #  x <- as.numeric(levels(x))[x]
    #}
    ##Function 2
    isTxt <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0

    ##Function 3
    ####function to extract dirInfo
    dirInfo <- function(x){
        str_replace_all(x, pattern = "\\\\", replacement="/") %>%
            str_split( "/") %>%
            unlist() %>%
            as.data.frame()  %>%
            mutate(level = row_number())
    }

    ##Function 4
    ########Processing Image
    procImage <- function(x){
        tmpImage <- readImage(x)#file_list.name
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
        write.table(singleData, paste0(x,"_",sdate,"_ImageData.txt"),
                    sep="\t",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
        message(sprintf("%s_%s_ImageData.txt",x,sdate))

        if (Binary == TRUE){
            dmrmabw <- 1*(bnrySeg > 0)# %>% # binary image
            file.name <- paste0(x,"_",sdate,  "_Binary.", imgType)
            writeImage(dmrmabw, file.name,type = imgType,quality = 100)
            message(file.name)

        }
        nDmrmadm <- nDm*(bnrySeg > 0)# %>% # distance map

        file.name <- paste0(x,"_",sdate,  "_DistMap.", imgType)
        writeImage(nDmrmadm, file.name,type = imgType,quality = 100)
        message(file.name)
    }



    ##Function 5
    ######Selecting the Directory#######
    procFiles <- function() {  ##recursive function
        filePath <- choose.files(caption = "Select any file to set the directory", multi = FALSE)

        if (length(filePath) == 0) return(message("\nFile not selected.\n"))

        # extract directory info from the selected file
        dirInfo1 <- dirInfo(filePath)
        colnames(dirInfo1) <- c("Dir_Name", "Level")
        setwd(str_c(dirInfo1$Dir_Name[seq_len(nrow(dirInfo1) - 1)], collapse = "/"))

        lapply(list.files()[isTxt(list.files())], procImage)

        procFiles()
    }

    ##Excecute
    procFiles()

}












