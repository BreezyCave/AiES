#' @title creating distance map and binary image from tiff image file
#' @description \code{Axdistmap} create images of distance map and binary image (option) from tiff image files
#' @import stringr
#' @import png
#' @import colorspace
#' @import EBImage
#' @import dplyr
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom stats na.omit
#' @param Binary TRUE: exporting binary image files
#' @param All_Features TRUE: exporting data of all features
#' @return return the image of distancemap and data of features
#' @export
#' @examples
#' # Axdistmap(Binary = FALSE, All_Features = FALSE)

Axdistmap <- function( Binary = FALSE, All_Features = FALSE){

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

  ####Information of Date####
  sdate <- Sys.Date()

  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }



  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0

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
      rma <- 30
      rmNum <- which(sdat$s.area <= rma)
      dml <- rmObjects(dml, rmNum)
      sdat <- as.data.frame(computeFeatures.shape(dml)) #option with
      ####this elimination option is not included in the original program
      #colorMode(dml) <- Grayscale

      hdat <- as.data.frame(computeFeatures.haralick(dml,ddm))
      mdat <- as.data.frame(computeFeatures.moment(dml))

      data_sht <<- data.frame(cbind( sdat, hdat, mdat))


      invisible({rm(list=c("sdat", "hdat", "mdat"));gc();gc()})
      data_sht <- na.omit(data_sht)
      if (is.numeric(data_sht[,1]) == FALSE){
        data_sht <- sapply(data_sht,num_f)
      }

      #sdat0 <- as.data.frame(computeFeatures.shape(dml))#When performing this step, an error may occur due to the file size being too large
      ##sdat <- as.data.frame(computeFeatures.shape(dmrmal))#When performing this step, an error may occur due to the file size being too large

      if (Binary == TRUE){
        dmrmabw <- 1*(dml > 0)# binary image
        file.name <- paste0(sdate, file_list.name, "_Binary.png")
        writePNG(dmrmabw,file.name)
      }
      ddmrmadm <- ddm*(dml > 0)# distance map
      file.name <- paste0(sdate, file_list.name, "_DistMap.png")
      writePNG(ddmrmadm,file.name)
      dmrmal <- bwlabel(ddmrmadm)
      colorMode(dmrmal) <- Grayscale

      cols = c('black', sample(heat_hcl(max(dmrmal))))
      dHeat = Image(cols[1+dmrmal], dim=dim(dmrmal))
      file.name <- paste0(sdate, file_list.name, "_Color.png")
      writePNG(dHeat,file.name)
      ##############################

      ##sdat <- as.data.frame(computeFeatures.shape(dmrmal))#When performing this step, an error may occur due to the file size being too large


      if (exists('data_sh') == FALSE) {
      } else {
        data_sh <- rbind(data_sh, data.frame(cbind("File" = file_list.name,  "Object size percentage" = sum(data_sht$s.area)/(dim(ddmrmadm)[1]*dim(ddmrmadm)[2])*100)))
      }
  }
      write.table(data_sh, paste0(sdate,"_Area.txt"),sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
      return(data_sh)
}
