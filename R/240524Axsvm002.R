#' @title creating a classifer based on SVM-based machine learning
#' @description \code{Axsvm} create a classifer
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
#' # Axsvm(Image = TRUE, nCst = 3, nGmm = 0.1, nCrss=5)


Axsvm <- function(nCst = 3, nGmm = 0.1, nCrss=5){

######################Codes for test run
#  library("dplyr")
#  library("e1071")
#  library("stringr")
#  library("ggpubr")
#  Image = TRUE
#  nCst = 3
#  nGmm = 0.1
#  nCrss=5
#######################################


  ##############Definition of variable######################

  txt_Group <- c("Degenerate","Intact")
  Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")
  sdate <- Sys.Date()

  ###############Definition of function######################
  is.txt <- function(x) regexpr('\\.txt$', x) > 0
  ####function to extract dir_info
  dir_info <- function(x){
    str_replace_all(x, pattern = "\\\\", replacement="/") %>%
      str_split( "/") %>%
      unlist() %>%
      as.data.frame()  %>%
      mutate(level = row_number())
  }
  ################################


  ######Selecting the Directory1#######
  # file full path
    file_path1 <- choose.files(caption = "Select a .txt file to set the directory (degenerate condition) ",
                               multi=FALSE)

    if (length(file_path1)==0) {
      cat("File not selected.\n")
      return(NULL)##Error countermeasure code
    }else{
    # extract directory info from the selected file
    dir_info1 <- dir_info(file_path1)
    colnames(dir_info1) <- c("Dir_Name", "Level")
    Mainf[1] <- str_c(dir_info1$Dir_Name[1:(nrow(dir_info1)-1)], collapse = "/")
    #####################################
    }

  ######Selecting the Directory2#######
  # file full path
  file_path2 <- choose.files(caption = "Select a .txt file to set the directory (intact condition) ",
                             multi=FALSE)

  if (length(file_path2)==0) {
    cat("File not selected.\n")
    return(NULL)##Error countermeasure code
  }else{
    # extract directory info from the selected file
  dir_info2 <- dir_info(file_path2)
  colnames(dir_info2) <- c("Dir_Name", "Level")
  Mainf[2] <- str_c(dir_info2$Dir_Name[1:(nrow(dir_info2)-1)], collapse = "/")
  #####################################
  }


###############Create training date####

  tryCatch({
    if(exists('Data_sh') == TRUE) invisible({rm(Data_sh);gc();gc()})
    }condition = function(c){
      cat("Object 'Data_sh' deleted\n")
    })

  for (i in 1:2) {
      setwd(Mainf[i])
    for (file_list.name in list.files()[is.txt(list.files())]){

      Temp_data <- read.table(file_list.name, header=T, sep="\t")

      if (exists('Data_sh') == FALSE) {
        Data_sh <- data.frame(cbind("Group" = as.factor(txt_Group[i]), Temp_data),
                              stringsAsFactors = TRUE)
      } else {
        Data_sh <- rbind(Data_sh, data.frame(cbind("Group" = as.factor(txt_Group[i]), Temp_data),
                                             stringsAsFactors = TRUE))
      }
    }
  }
########################################

##################delete unnecssary columns
  data_svm <- Data_sh[ which(colnames(Data_sh) %in% Ftr.tmp)]
###########################################


  ##########################################
  #SVM
    SVM_model <- svm(
      Group ~ ., data = data_svm,
      method = "C-classification",
      probability = TRUE,
      kernel = "radial",
      gamma = nGmm,
      cost = nCst,
      cross = nCrss, #K-fold cross-validation
    )
  SVM_model

  ##########################export SVM model file
  ######Selecting the Directory1#######
  # file full path
  file_path4 <- choose.files(default = paste0(sdate,"_AxClassifer.svm") ,caption = "Save a SVM model", multi = FALSE)

  if (length(file_path4)==0) {
    cat("Canceled \n")
    return(NULL)##Error countermeasure code
  }else{
    save(SVM_model, file= file_path4)
  }

  ######Export extracted data file#######
  ######Selecting the Directory1#######
  # file full path
  file_path3 <- choose.files(default = paste0(sdate,"_Extracted_data_for_ML.txt") ,caption = "Save a extracted data file", multi = FALSE)
  dir_info3 <- dir_info(file_path3)
  colnames(dir_info3) <- c("Dir_Name", "Level")

  if (length(file_path3)==0) {
    cat("Canceled \n")
    return(NULL)##Error countermeasure code
  }else{
    write.table(data_svm, file_path3,
                sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    print(dir_info3$Dir_Name[nrow(dir_info3)])
  }





}###end of function








