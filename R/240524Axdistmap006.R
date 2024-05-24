#' @title creating distance map and binary image from tiff image file
#' @description \code{Axdistmap} create images of distance map and binary image (option) from tiff image files
#' @import stringr
#' @import colorspace
#' @import dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom stats na.omit
#' @param Sub_Back 30: subtract background objects (default 30 pixels)
#' @param Binary TRUE: exporting binary image files
#' @param All_Features TRUE: exporting data of all features
#' @param Type png, jpg, tiff
#' @return return the image of distancemap and data of features
#' @export
#' @examples
#' # Axdistmap(Sub_Bacl = 30, Binary = FALSE, All_Features = FALSE, Type = tiff)


Axdistmap <- function(Sub_Back = 30, Binary = FALSE, All_Features = FALSE, Type = "tiff"){

  #### #test code for debug
  #library("EBImage")
  #library("stringr")
  #library("colorspace")
  #library("EBImage")
  #library("dplyr")
  #Sub_Back <- 30
  #Binary = "TRUE"
  #All_Features = "TRUE"
  #Type = "tiff"
  #########################
  ###############Definition of variable######################

  sdate <- Sys.Date()

  ###############Definition of function######################
  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }
  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0
  ####function to extract dir_info
  dir_info <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }


  ######Selecting the Directory#######
  file_path <- choose.files(caption = "Select any file to set the directory",
                            multi=FALSE)


  if (length(file_path)==0) {
    cat("File not selected.\n")
    return()##Error countermeasure code
  }
  while(length(file_path)!=0){
    # extract directory info from the selected file
    dir_info <- dir_info(file_path)
    colnames(dir_info) <- c("Dir_Name", "Level")
    setwd(str_c(dir_info$Dir_Name[1:(nrow(dir_info)-1)], collapse = "/"))


    for (file_list.name in list.files()[is.tif(list.files())]){

      test <- readImage(file_list.name)
      test <- resize(test, w = 900)#15Jun11:696 22Apr15:900
      test_s <- (1/(1+(0.5/test[,,1])^5))
      mf1 <- medianFilter(test_s,1)
      x_mf1 <- (log1p(mf1)/log1p(max(mf1)))
      x_mf1_log <- (log1p(x_mf1)/log1p(max(x_mf1)))
      x_thr <- thresh(x_mf1_log, 8, 8, -0.04)
      kern <- makeBrush(1, shape="box")
      x_thr_bw <- 1 - opening(closing(x_thr, kern), kern)
      dm <- distmap(x_thr_bw)
      ddm <- normalize(dm)
      dml <- bwlabel(ddm)
      sdat <- as.data.frame(computeFeatures.shape(dml))
      #######Threshold for eliminating small objects that cannot be classified by image processing
      rma <- Sub_Back
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

      if(All_Features == FALSE){
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
      write.table(data_sht, paste0(file_list.name,sdate,"_ImageData.txt"),
                  sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
      print(paste0(file_list.name,sdate,"_ImageData.txt"))

      if (Binary == TRUE){
        dmrmabw <- 1*(dml > 0)# %>% # binary image
        if(Type == c("png")){
          file.name <- paste0( file_list.name,sdate, "_Binary.png")
          writeImage(dmrmabw, file.name,type = "png",quality = 100)
          print(file.name)
        }else if(Type == c("jpg")){
          file.name <- paste0(file_list.name,sdate,  "_Binary.jpg")
          writeImage(dmrmabw, file.name,type = "jpg",quality = 100)
          print(file.name)
        }else{
          file.name <- paste0( file_list.name,sdate, "_Binary.tiff")
          writeImage(dmrmabw, file.name,type = "tiff",quality = 100)
          print(file.name)
        }
      }
      ddmrmadm <- ddm*(dml > 0)# %>% # distance map
      if(Type == c("png")){
        file.name <- paste0(file_list.name,sdate,  "_DistMap.png")
        writeImage(ddmrmadm, file.name,type = "png",quality = 100)
        print(file.name)
      }else if(Type == c("jpg")){
        file.name <- paste0( file_list.name,sdate, "_DistMap.jpg")
        writeImage(ddmrmadm, file.name,type = "jpg",quality = 100)
        print(file.name)
      }else{
        file.name <- paste0( file_list.name,sdate, "_DistMap.tiff")
        writeImage(ddmrmadm, file.name,type = "tiff",quality = 100)
        print(file.name)
      }

    }
    file_path <- choose.files(caption = "Select any file to set the directory",
                                multi=FALSE)
  }
}
