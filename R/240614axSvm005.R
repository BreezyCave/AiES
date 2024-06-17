# Copyright 2024 Your Company Name
# BSD 3-Clause License (see LICENSE file)
#' @title creating a classifer based on SVM-based machine learning
#' @description \code{axSvm} create a classifer
#' @import stringr
#' @import dplyr
#' @importFrom e1071 svm
#' @importFrom ggpubr mutate
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom utils read.table
#' @param nCst 3: cost parameter for libsvm
#' @param nGmm 0.1: gamma parameter for libsvm
#' @param nCrss 5: number of K-fold cross-validation
#' @return return the svm model and test data
#' @export
#' @examples
#' # axSvm(Image = TRUE, nCst = 3, nGmm = 0.1, nCrss=5)


axSvm <- function(nCst = 3, nGmm = 0.1, nCrss=5){

######################Codes for test run
  ##library("dplyr")
  ##library("e1071")
  ##library("stringr")
  ##library("ggpubr")
  ##Image = TRUE
  ##nCst = 3
  ##nGmm = 0.1
  ##nCrss=5
#######################################


  ##############Definition of variable######################

  txtGroup <- c("Degenerate","Intact")
  ftrTmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  sdate <- Sys.Date()

  ###############Definition of function######################
  ##Function 1
  createEmptyDF <- function( nrow, ncol, colnames = c() ){
    if( missing( ncol ) && length( colnames ) > 0 ){
      ncol <- length( colnames )
    }
    data.frame( matrix( vector(), nrow, ncol, dimnames = list( c(), colnames ) ) )
  }

  ##Function2
  isTxt <- function(x) regexpr('\\.txt$', x) > 0
  ####function to extract dirInfo
  dirInfo <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }
  ################################


  ######Selecting the Directory1#######
  # file full path
    filePath1 <- choose.files(caption = "Select a .txt file to set the directory (degenerate condition) ",
                               multi=FALSE)

    if (length(filePath1)==0) {
      return(message("File not selected.\n"))##Error countermeasure code
    }else{
    # extract directory info from the selected file
    dirInfo1 <- dirInfo(filePath1)
    colnames(dirInfo1) <- c("Dir_Name", "Level")
    mainFldr[1] <- str_c(dirInfo1$Dir_Name[seq_len(nrow(dirInfo1)-1)], collapse = "/")
    #####################################
    }

  ######Selecting the Directory2#######
  # file full path
  filePath2 <- choose.files(caption = "Select a .txt file to set the directory (intact condition) ",
                             multi=FALSE)

  if (length(filePath2)==0) {
    return(message("File not selected.\n"))##Error countermeasure code
  }else{
    # extract directory info from the selected file
  dirInfo2 <- dirInfo(filePath2)
  colnames(dirInfo2) <- c("Dir_Name", "Level")
  mainFldr[2] <- str_c(dirInfo2$Dir_Name[seq_len(nrow(dirInfo2)-1)], collapse = "/")
  #####################################
  }


###############Create training date####

  ##tryCatch({
  ##  if(exists('trainingData') == TRUE) invisible({rm(trainingData);gc();gc()})
  ##}, condition = function(c){
  ##  warning("Object 'trainingData' deleted\n")
  ##})

  ##

  ##Function
  procFiles <- function(x,y){
    TempData <- read.table(x, header=TRUE, sep="\t")

    if (exists('trainingData') == FALSE) {
      trainingData <- data.frame(cbind("Group" = as.factor(txtGroup[y]), TempData),
                            stringsAsFactors = TRUE)
    } else {
      trainingData <- rbind(trainingData, data.frame(cbind("Group" = as.factor(txtGroup[y]), TempData),
                                           stringsAsFactors = TRUE))
    }
    return(trainingData)
  }


  ##Function
  procSvm <- function(x){
    setwd(mainFldr[x])
    dataList <- lapply(list.files()[isTxt(list.files())], procFiles, y=x)
    dataList <- do.call(rbind, dataList)

    ##################delete unnecssary columns
    dataList <- dataList[ which(colnames(dataList) %in% ftrTmp)]
    ###########################################

    return(dataList)

  }

  ####################Execution
    dataSvm <- lapply(seq(txtGroup), FUN = procSvm)
    dataSvm <- do.call(rbind, dataSvm)

    ##########################################
    #SVM
    svmModel <- svm(
      Group ~ ., data = dataSvm,
      method = "C-classification",
      probability = TRUE,
      kernel = "radial",
      gamma = nGmm,
      cost = nCst,
      cross = nCrss, #K-fold cross-validation
    )
    svmModel


    ######Export extracted data file#######
    ######Selecting the Directory1#######
    # file full path
    filePath3 <- choose.files(default = paste0(sdate,"_Extracted_data_for_ML.txt") ,caption = "Save a extracted data file", multi = FALSE)
    dirInfo3 <- dirInfo(filePath3)
    colnames(dirInfo3) <- c("Dir_Name", "Level")

    if (length(filePath3)==0) {
      message("Canceled \n")
      return(NULL)##Error countermeasure code
    }else{
      write.table(data_svm, filePath3,
                  sep="\t",row.names=FALSE, quote=FALSE, col.names=TRUE, append=FALSE)
      message(dirInfo3$Dir_Name[nrow(dirInfo3)])
    }

    ##########################export SVM model file
    ######Selecting the Directory1#######
    # file full path
    filePath4 <- choose.files(default = paste0(sdate,"_AxClassifer.svm") ,caption = "Save a SVM model", multi = FALSE)
    dirInfo4 <- dirInfo(filePath4)
    colnames(dirInfo4) <- c("Dir_Name", "Level")

    if (length(filePath4)==0) {
      message("Canceled \n")
      return(NULL)##Error countermeasure code
    }else{
      save(svmModel, file= filePath4)
      message(dirInfo4$Dir_Name[nrow(dirInfo4)])
    }



##  for (i in 1:2) {
##      setwd(mainFldr[i])
    ##    for (file_list.name in list.files()[isTxt(list.files())]){
    ##
    ##TempData <- read.table(file_list.name, header=TRUE, sep="\t")
    ##
    ##if (exists('trainingData') == FALSE) {
    ##  trainingData <- data.frame(cbind("Group" = as.factor(txtGroup[i]), TempData),
    ##                        stringsAsFactors = TRUE)
    ##} else {
    ##  trainingData <- rbind(trainingData, data.frame(cbind("Group" = as.factor(txtGroup[i]), TempData),
    ##                                       stringsAsFactors = TRUE))
    ##}
    ##}
##  }
########################################



}###end of function








