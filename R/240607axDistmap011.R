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
  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }
  ##Function 2
  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0

  ##Function 3
  ####function to extract dir_info
  dir_info <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }

  ##Function 4
  ########Processing Image
  procImage <- function(x){
    test <- readImage(x)#file_list.name
    test <- resize(test, w = 900)#15Jun11:696 22Apr15:900
    test_s <- (1/(1+(0.5/test[,,1])^5))
    mf1 <- medianFilter(test_s,1)
    x_mf1 <- (log1p(mf1)/log1p(max(mf1)))
    x_mf1_log <- (log1p(x_mf1)/log1p(max(x_mf1)))
    x_thr <- 1 - thresh(x_mf1_log, 8, 8, -0.04)
    kern <- makeBrush(1, shape="box")
    x_thr_bw <- opening(closing(x_thr, kern), kern)
    dm <- distmap(x_thr_bw)
    ddm <- normalize(dm)
    dml <- bwlabel(ddm)
    sdat <- as.data.frame(computeFeatures.shape(dml))
  #######Threshold for eliminating small objects that cannot be classified by image processing
    rma <- subBack
    rmNum <- which(sdat$s.area <= rma)
    dml <- rmObjects(dml, rmNum)
    sdat <- as.data.frame(computeFeatures.shape(dml)) #option with
    hdat <- as.data.frame(computeFeatures.haralick(dml,ddm))
    mdat <- as.data.frame(computeFeatures.moment(dml))
    sing_data <- data.frame(cbind(sdat, hdat, mdat))
    invisible({rm(list=c("sdat", "hdat", "mdat"));gc();gc()})
    data_sht <- na.omit(sing_data)
    if (is.numeric(data_sht[,1]) == FALSE){
      data_sht <- sapply(data_sht,num_f)
    }

    if(allFeatures == FALSE){
      data_sht <- data.frame(s.area = data_sht$s.area,
                            m.eccentricity = data_sht$m.eccentricity,
                            s.radius.sd = data_sht$s.radius.sd,
                            h.sva.s2 = data_sht$h.sva.s2,
                            h.idm.s1 = data_sht$h.idm.s1,
                            h.sen.s1 = data_sht$h.sen.s1,
                            m.majoraxis = data_sht$m.majoraxis,
                            stringsAsFactors = TRUE)
    }else {
      data_sht <- data.frame(cbind(data_sht,Cir = data_sht$s.area*pi*4/data_sht$s.perimeter^2) , stringsAsFactors = TRUE)
    }

  #####export data file
    write.table(data_sht, paste0(x,"_",sdate,"_ImageData.txt"),
                sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    print(paste0(x,sdate,"_ImageData.txt"))

    if (Binary == TRUE){
      dmrmabw <- 1*(dml > 0)# %>% # binary image
      file.name <- paste0(x,"_",sdate,  "_Binary.", imgType)
      writeImage(dmrmabw, file.name,type = imgType,quality = 100)
      print(file.name)

    }
    ddmrmadm <- ddm*(dml > 0)# %>% # distance map

    file.name <- paste0(x,"_",sdate,  "_DistMap.", imgType)
    writeImage(ddmrmadm, file.name,type = imgType,quality = 100)
    print(file.name)
  }



  ##Function 5
  ######Selecting the Directory#######
  process_files <- function() {  ##recursive function
    file_path <- choose.files(caption = "Select any file to set the directory", multi = FALSE)

    if (length(file_path) == 0) return(cat("\nFile not selected.\n"))

    # extract directory info from the selected file
    dir_info1 <- dir_info(file_path)
    colnames(dir_info1) <- c("Dir_Name", "Level")
    setwd(str_c(dir_info1$Dir_Name[1:(nrow(dir_info1) - 1)], collapse = "/"))

    lapply(list.files()[is.tif(list.files())], procImage)

    process_files()
  }

  ##Excecute
  process_files()

}












