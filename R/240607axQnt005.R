# Copyright 2024 Your Company Name
# BSD 3-Clause License (see LICENSE file)
#' @title Axon integrity index quantification
#' @description \code{axQnt} quantificate axonal integrity index (and degeneration index) and export data file
#' @import stringr
#' @import dplyr
#' @importFrom EBImage readImage resize medianFilter thresh makeBrush opening closing distmap bwlabel rmObjects normalize computeFeatures.shape computeFeatures.haralick computeFeatures.moment
#' @importFrom stats predict
#' @importFrom ggpubr mutate
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom utils read.table
#' @importFrom stats na.omit
#' @param imprtImg TRUE: import image file (.tiff) or data file (.txt)
#' @param subBack 30: subtract background objects (default 30 pixels)
#' @param expSip TRUE: export single image prediction data file (.csv)
#' @return return the Summary csv data and Image data (option)
#' @export
#' @examples
#' # axQnt(imprtImg = TRUE, subBack = 30, expSip = TRUE)


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

  txt_Group <- c("Degenerate","Intact")
  Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  Ftr.svm <- c("m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  Rlt.name <- c("FileName","AxII","DegI")
  sdate <- Sys.Date()

  ###############Definition of function######################
  ##Function 1
  CreateEmptyDF = function( nrow, ncol, colnames = c() ){
    if( missing( ncol ) && length( colnames ) > 0 ){
      ncol = length( colnames )
    }
    data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
  }

  ##Function 2
  is.txt <- function(x) regexpr('\\.txt$', x) > 0
  ##Function 3
  is.tif <- function(x) regexpr('\\.tif$', x) + regexpr('\\.tiff$', x)> 0
  ##Function 4
  num_f <- function(x){
    x <- as.numeric(levels(x))[x]
  }

  ##Function 5
  ####function to extract dir_info
  dir_info <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }

  ##Function 6
  ########Processing Image
  procImage <- function(x){
    test <- readImage(x)#file_list.name
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
    rma <- subBack
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
    write.table(data_sht, paste0(x,sdate,"_ImageData.txt"),
                sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    print(paste0(x,"_",sdate,"_ImageData.txt"))
  }

  ##Function 7
  ########Processing Text data
  procTxt <- function(x, y){
    Data_shq <- read.table(x, header=T, sep="\t")
    if(! FALSE %in% (Ftr.svm %in% colnames(Data_shq))){
      Pred <- svm_model_loaded %>%
      get() %>%
      predict(Data_shq, type="class", probability = FALSE) %>%
      as.character()
      if(expSip == TRUE){
        write.table(data.frame(Pred = Pred,Data_shq),
                  paste0(x, "_", sdate,"_SIP.csv"),
                  sep=",",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
        print(paste0(x, "_", sdate,"_SIP.csv"))
      }
      AxII <- sum((Pred %in% txt_Group[2])*Data_shq$s.area)/(sum((Pred %in% txt_Group[2])*Data_shq$s.area)+sum((Pred %in% txt_Group[1])*Data_shq$s.area))
      DegI <- sum((Pred %in% txt_Group[1])*Data_shq$s.area)/(sum((Pred %in% txt_Group[2])*Data_shq$s.area)+sum((Pred %in% txt_Group[1])*Data_shq$s.area))
      y[nrow(y)+1,] <- c(x, AxII, DegI)
      return(y)
    }else{
      cat(x, ": This file does not contain the required data \n")
    }
  }

  ##Function 8
  ######Selecting the Directory#######
  process_files <- function() {  ##recursive function
    # file full path
    ifelse(imprtImg == TRUE,
           file_path2 <- choose.files(caption = "Select any tiff file to set the directory",
                                      multi=FALSE),
           file_path2 <- choose.files(caption = "Select any text file to set the directory",
                                      multi=FALSE))
    if(length(file_path2)==0) return(cat("File not selected.\n"))

    # extract directory info from the selected file
    dir_info2 <- dir_info(file_path2)
    colnames(dir_info2) <- c("Dir_Name", "Level")
    setwd(str_c(dir_info2$Dir_Name[1:(nrow(dir_info2)-1)], collapse = "/"))


    if(imprtImg == TRUE){
      ################Calculate axon integrity index from .tiff data#####
      lapply(list.files()[is.tif(list.files())], procImage)
    }

    ##
    ################Calculate axon integrity index from .txt data#####
    Rlt.Summary <- CreateEmptyDF(0, colnames = Rlt.name)
    Rlt.Summary <- lapply(list.files()[is.txt(list.files())], procTxt, y=Rlt.Summary)
    Rlt.Summary <- do.call(rbind, Rlt.Summary)

    ##################Export the Summary data (csv)##############################
    write.table(Rlt.Summary, paste0(dir_info2$Dir_Name[(nrow(dir_info2)-1)], "_", sdate,"_Summary.csv"),
                sep=",",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    print(paste0(dir_info2$Dir_Name[(nrow(dir_info2)-1)], "_", sdate,"_Summary.csv"))


    process_files()##recursive function
  }

  ################################
  ##Excecute
  #############import svm model #######
  ######Selecting the Directory1#######
  # file full path
  file_path1 <- choose.files(caption = "Select a SVM model",
                             multi=FALSE)
  if (length(file_path1)==0) {
    return(cat("File not selected.\n"))##Error countermeasure code
  }else{
    tryCatch({
      dir_info1 <- dir_info(file_path1)
      colnames(dir_info1) <- c("Dir_Name", "Level")
      svm_model_loaded = load(file_path1)
    }, error = function(e){
      cat(" Error in load file '",
          dir_info1$Dir_Name[nrow(dir_info1)],"'\n This file dose not have SVM model. \n")
    })
  }

  #####################################
  process_files()



}#####end of function








