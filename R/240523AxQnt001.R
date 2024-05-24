#' @title Axon integrity index quantification
#' @description \code{AxQnt} quantificate axonal integrity index (and degeneration index) and export data file
#' @import stringr
#' @import dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment
#' @importFrom stats predict
#' @importFrom ggpubr mutate
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom utils read.table
#' @importFrom stats na.omit
#' @param imprt_IMG TRUE: import image file (.tiff) or data file (.txt)
#' @param Sub_Back 30: subtract background objects (default 30 pixels)
#' @param exp_SIP TRUE: export single image prediction data file (.csv)
#' @return return the Summary csv data and Image data (option)
#' @export
#' @examples
#' # AxQnt(imprt_IMG = TRUE, Sub_Back = 30, exp_SIP = TRUE)


AxQnt <- function(imprt_IMG = TRUE, Sub_Back = 30, exp_SIP = TRUE){

  ######################Codes for test run
  #library("dplyr")
  #library("stats")
  #library("stringr")
  #library("ggpubr")
  #library("EBImage")
  #library("stringr")
  #library("colorspace")
  #library("EBImage")
  #library("dplyr")
  #Sub_Back <- 30
  #imprt_IMG = TRUE
  #exp_SIP = TRUE
  ###############Definition of variable######################

  txt_Group <- c("Degenerate","Intact")
  Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  Ftr.svm <- c("m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  Rlt.name <- c("FileName","AxII","DegI")

  ###############Definition of function######################
  CreateEmptyDF = function( nrow, ncol, colnames = c() ){
    if( missing( ncol ) && length( colnames ) > 0 ){
      ncol = length( colnames )
    }
    data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
  }



########################################
  sdate <- Sys.Date()

  is.txt <- function(x) regexpr('\\.txt$', x) > 0
  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0
  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }

  ####function to extract dir_info
  dir_info <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }
  ################################


  #############import svm model #######
  ######Selecting the Directory1#######
  # file full path
  file_path1 <- choose.files(caption = "Select a SVM model",
                               multi=FALSE)
  dir_info1 <- dir_info(file_path1)
  colnames(dir_info1) <- c("Dir_Name", "Level")
  if (length(file_path1)==0) {
    cat("File not selected.\n")
    return()##Error countermeasure code
  }else{
    tryCatch({
      svm_model_loaded = load(file_path1)
    }, error = function(e){
      cat(" Error in load file '",
          dir_info1$Dir_Name[nrow(dir_info1)],"'\n This file dose not have SVM model. \n")
    })
  }

  #####################################

  ######Selecting the Directory2####################################
  # file full path
  if(imprt_IMG == TRUE){
  file_path2 <- choose.files(caption = "Select any tiff file to set the directory",
                             multi=FALSE)
  }else{
  file_path2 <- choose.files(caption = "Select any text file to set the directory",
                             multi=FALSE)
  }
  if (length(file_path2)==0) {
    cat("File not selected.\n")
    return()##Error countermeasure code
  }
  while(length(file_path2)!=0){
    # extract directory info from the selected file
    dir_info2 <- dir_info(file_path2)
    colnames(dir_info2) <- c("Dir_Name", "Level")
    setwd(str_c(dir_info2$Dir_Name[1:(nrow(dir_info2)-1)], collapse = "/"))

    Rlt.Summary <- CreateEmptyDF(0, colnames = Rlt.name)


    if(imprt_IMG == TRUE){


      ################Calculate axon integrity index from .tiff data#####
      ###############################################
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
        sdat <- as.data.frame(computeFeatures.shape(dml))
        hdat <- as.data.frame(computeFeatures.haralick(dml,ddm))
        mdat <- as.data.frame(computeFeatures.moment(dml))

        sing_data <- data.frame(cbind(sdat, hdat, mdat))
        invisible({rm(list=c("sdat", "hdat", "mdat"));gc();gc()})
        data_sht <- na.omit(sing_data)
        if (is.numeric(data_sht[,1]) == FALSE){
          data_sht <- sapply(data_sht,num_f)
        }

        data_sht <- data.frame(s.area = data_sht$s.area,
                               m.eccentricity = data_sht$m.eccentricity,
                               s.radius.sd = data_sht$s.radius.sd,
                               h.sva.s2 = data_sht$h.sva.s2,
                               h.idm.s1 = data_sht$h.idm.s1,
                               h.sen.s1 = data_sht$h.sen.s1,
                               m.majoraxis = data_sht$m.majoraxis,
                               stringsAsFactors = TRUE)

        #####export data file
        write.table(data_sht, paste0(file_list.name,sdate,"_ImageData.txt"),
                    sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
      }
    }

################Calculate axon integrity index from .txt data#####
###############################################import txt data######
    for (file_list.name in list.files()[is.txt(list.files())]){
      Data_sh <- read.table(file_list.name, header=T, sep="\t")
      if(! FALSE %in% (Ftr.svm %in% colnames(Data_sh))){
        Pred <- svm_model_loaded %>%
          predict(Data_sh, type="class", probability = FALSE) %>%
          as.character()
        if(exp_SIP == TRUE){
          write.table(data.frame(Pred = Pred,Data_sh),
                    paste0(file_list.name, "_", sdate,"_SIP.csv"),
                    sep=",",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
        }
        AxII <- sum((Pred %in% txt_Group[2])*Data_sh$s.area)/(sum((Pred %in% txt_Group[2])*Data_sh$s.area)+sum((Pred %in% txt_Group[1])*Data_sh$s.area))
        DegI <- sum((Pred %in% txt_Group[1])*Data_sh$s.area)/(sum((Pred %in% txt_Group[2])*Data_sh$s.area)+sum((Pred %in% txt_Group[1])*Data_sh$s.area))
        Rlt.Summary[nrow(Rlt.Summary)+1,] <- c(file_list.name, AxII, DegI)
      }else{
        cat(file_list.name, ": This file does not contain the required data \n")
      }

    }


##################Export the Summary data (csv)##############################
    write.table(Rlt.Summary, paste0(dir_info2$Dir_Name[(nrow(dir_info2)-1)], "_", sdate,"_Summary.csv"),
                sep=",",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)


##################change the directory##############################
    if(imprt_IMG == TRUE){
      file_path2 <- choose.files(caption = "Select any tiff file to set the directory",
                                 multi=FALSE)
    }else{
      file_path2 <- choose.files(caption = "Select any text file to set the directory",
                                 multi=FALSE)
    }


    if (length(file_path2)==0) {
      cat("File not selected.\n")
      return()##Error countermeasure code
    }
  }

}#####end of function








