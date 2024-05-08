#' @title creating distance map and binary image from tiff image file
#' @description \code{Axdistmap} create images of distance map and binary image (option) from tiff image files
#' @import stringr
#' @import colorspace
#' @import dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment writeImage
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom stats na.omit
#' @param Sub_Back TRUE: subtract background under 30 pixels
#' @param Binary TRUE: exporting binary image files
#' @param All_Features TRUE: exporting data of all features
#' @param Type png, jpg, tiff
#' @return return the image of distancemap and data of features
#' @export
#' @examples
#' # Axdistmap(Sub_Bacl = TRUE, Binary = FALSE, All_Features = FALSE, Type = tiff)


Axdistmap <- function(Sub_Back = TRUE, Binary = FALSE, All_Features = FALSE, Type = "tiff"){

  #### #test code for debug
  #library("EBImage")
  #library("stringr")
  #library("colorspace")
  #library("EBImage")
  #library("dplyr")
  #Sub_Back = "FALSE"
  #Binary = "TRUE"
  #All_Features = "TRUE"
  #Type = "tiff"
  #########################

  sdate <- Sys.Date()
  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }
  ######Selecting the Directory#######
  # 代表的なファイルのファイルパス
  file_path <- choose.files()

  # ディレクトリ情報の抽出
  dir_info <- str_replace_all(file_path, pattern = "\\\\", replacement="/") %>%
    str_split( "/") %>%
    unlist() %>%
    as.data.frame()  %>%
    mutate(level = row_number())
  colnames(dir_info) <- c("Dir_Name", "Level")
  # 最後のファイル名を除いた部分を抽出
  Mainf <- str_c(dir_info$Dir_Name[1:(nrow(dir_info)-1)], collapse = "/")
  setwd(Mainf)#Low probability of error due to direct selecting a file through dialogue box
  ######Selecting the Directory#######


  #####################
  # ディレクトリ情報の抽出
  #file_info <- str_replace_all(file_list.name, pattern = "\\\\", replacement="/") %>%
  #str_split( "/") %>%
  #unlist() %>%
  #as.data.frame()  %>%
  #mutate(level = row_number())
  #colnames(file_info) <- c("Dir_Name", "Level")
  ###################

  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0

  ####Information of Date####


    for (file_list.name in list.files()[is.tif(list.files())]){
      #####################
      #test code for debug
      #sdate <- Sys.Date()

      #num_f <- function(x){
      #x <- as.numeric(levels(x))[x]
      #}

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
      if(Sub_Back == TRUE){
      rma <- 30
      rmNum <- which(sdat$s.area <= rma)
      dml <- rmObjects(dml, rmNum)
      sdat <- as.data.frame(computeFeatures.shape(dml)) #option with
      }
      ####this elimination option is not included in the original program
      #colorMode(dml) <- Grayscale
      ######Required ddm for caluculating haralick texture feature.
      ######However, the numbers of objects are different between them.


      hdat <- as.data.frame(computeFeatures.haralick(dml,ddm))
      mdat <- as.data.frame(computeFeatures.moment(dml))

      #
      #
      #sing_data <- data.frame(cbind("FileName" = file_list.name, sdat, hdat, mdat))
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
        data_sht <<- data.frame(cbind(data_sht,Cir = data_sht$s.area*pi*4/data_sht$s.perimeter^2) , stringsAsFactors = TRUE)

      }

      #if (nrow(sing_data)> NumPic) {
      #  sing_data <- sing_data[sample(nrow(sing_data), NumPic),]
      #}
      #if (exists('data_sh') == FALSE) {
      #  data_sh <<- sing_data
      #} else {
      #  data_sh <<- rbind(data_sh, sing_data)
      #}




      #####export data file
      write.table(data_sht, paste0(file_list.name,sdate,"_AllFeatures.txt"),
                  sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
      #return(data_sht)


      #sdat0 <- as.data.frame(computeFeatures.shape(dml))#When performing this step, an error may occur due to the file size being too large
      ##sdat <- as.data.frame(computeFeatures.shape(dmrmal))#When performing this step, an error may occur due to the file size being too large

      if (Binary == TRUE){
        dmrmabw <- 1*(dml > 0)# %>% # binary image
        #rotate( angle = -90) %>%
        #flop()
        if(Type == c("png")){
          file.name <- paste0( file_list.name,sdate, "_Binary.png")
          writeImage(dmrmabw, file.name,type = "png",quality = 100)
        }else if(Type == c("jpg")){
          file.name <- paste0(file_list.name,sdate,  "_Binary.jpg")
          writeImage(dmrmabw, file.name,type = "jpg",quality = 100)
        }else{
          file.name <- paste0( file_list.name,sdate, "_Binary.tiff")
          writeImage(dmrmabw, file.name,type = "tiff",quality = 100)
        }
      }
      ddmrmadm <- ddm*(dml > 0)# %>% # distance map
      #rotate( angle = -90) %>%
      #flop()
      if(Type == c("png")){
        file.name <- paste0(file_list.name,sdate,  "_DistMap.png")
        writeImage(ddmrmadm, file.name,type = "png",quality = 100)
      }else if(Type == c("jpg")){
        file.name <- paste0( file_list.name,sdate, "_DistMap.jpg")
        writeImage(ddmrmadm, file.name,type = "jpg",quality = 100)
      }else{
        file.name <- paste0( file_list.name,sdate, "_DistMap.tiff")
        writeImage(ddmrmadm, file.name,type = "tiff",quality = 100)
      }


      #dmrmal <- bwlabel(ddmrmadm)
      #colorMode(dmrmal) <- Grayscale

      #cols = c('black', sample(heat_hcl(max(dmrmal))))
      #dHeat = Image(cols[1+dmrmal], dim=dim(dmrmal))
      #file.name <- paste0(sdate, file_list.name, "_Color.png")
      #writePNG(dHeat,file.name)
      ##############################

      ##sdat <- as.data.frame(computeFeatures.shape(dmrmal))#When performing this step, an error may occur due to the file size being too large


      #if (exists('data_sh') == FALSE) {
      #} else {
      #  data_sh <- rbind(data_sh, data.frame(cbind("File" = file_list.name,  "Object size percentage" = sum(data_sht$s.area)/(dim(ddmrmadm)[1]*dim(ddmrmadm)[2])*100)))
      #}
  }
      #write.table(data_sh, paste0(sdate,"_Area.txt"),sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
      #return(data_sh)
}
