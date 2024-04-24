#' @title creating distance map and binary image from tiff image file
#' @description \code{Axdistmap} create images of distance map and binary image (option) from tiff image files
#' @import EBImage
#' @importFrom png writePNG
#' @importFrom colorspace heat_hcl
#' @importFrom utils write.table
#' @param x full path or the directory name containing tiff image files of phase contrast
#' @param CD TRUE: setting the directory "x" and "y" under the current working directory
#' @param Binary TRUE: exporting binary image files
#' @return return the list of file names and % area of axons
#' @export
#' @examples
#' # Axdistmap("C:/Users/R", Binary = TRUE) or Axdistmap("R", CD = TRUE, Binary = TRUE)


Axdistmap <- function(x, CD = FALSE, Binary = FALSE){
    Mainf <- x
    if (CD == TRUE){
      cdd <- getwd()
      Mainf <- paste(cdd, "/", x, sep="")
    }

    sdate <- Sys.Date()
    setwd(Mainf)#error 対処必要

    is.tif <- function(x) regexpr('\\.tif$', x) > 0
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
        colorMode(dml) <- Grayscale
        sdat0 <- as.data.frame(computeFeatures.shape(dml))#When performing this step, an error may occur due to the file size being too large
        rma <- 30 #Threshold for eliminating small objects that cannot be classified by image processing
        rmNum <- which(sdat0$s.area <= rma)
        dmlrma <- rmObjects(dml, rmNum)

        if (Binary == TRUE){
          dmrmabw <- 1*(dmlrma > 0)# binary image
          file.name <- paste0(sdate, file_list.name, "_Binary.png")
          writePNG(dmrmabw,file.name)
        }
        ddmrmadm <- ddm*(dmlrma > 0)# distance map
        file.name <- paste0(sdate, file_list.name, "_DistMap.png")
        writePNG(ddmrmadm,file.name)

        dmrmal <- bwlabel(ddmrmadm)
        colorMode(dmrmal) <- Grayscale

        cols = c('black', sample(heat_hcl(max(dmrmal))))
        dHeat = Image(cols[1+dmrmal], dim=dim(dmrmal))
        file.name <- paste0(sdate, file_list.name, "_Color.png")
        writePNG(dHeat,file.name)


        sdat <- as.data.frame(computeFeatures.shape(dmrmal))#When performing this step, an error may occur due to the file size being too large

        if (exists('data_sh') == FALSE) {
          data_sh <- data.frame(cbind("File" = file_list.name, "Object size percentage" = sum(sdat$s.area)/(dim(ddmrmadm)[1]*dim(ddmrmadm)[2])*100))
        } else {
          data_sh <- rbind(data_sh, data.frame(cbind("File" = file_list.name,  "Object size percentage" = sum(sdat$s.area)/(dim(ddmrmadm)[1]*dim(ddmrmadm)[2])*100)))
        }
  }
    write.table(data_sh, paste0(sdate,"_Area.txt"),sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    return(data_sh)
}
