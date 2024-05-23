#' @title creating a classifer based on SVM-based machine learning
#' @description \code{Axsvm} create a classifer
#' @import stringr
#' @import dplyr
#' @importFrom e1071 svm
#' @importFrom ggpubr mutate
#' @importFrom utils write.table
#' @importFrom utils choose.files
#' @importFrom utils read.table
#' @param Image TRUE: subtract background objects (under 30 pixels)
#' @param nCst 3: cost parameter for libsvm
#' @param nGmm 0.1: gamma parameter for libsvm
#' @param nCrss 5: number of K-fold cross-validation
#' @return return the svm model and test data
#' @export
#' @examples
#' # Axsvm(Image = TRUE, nCst = 3, nGmm = 0.1, nCrss=5)


Axsvm <- function(Image = TRUE, nCst = 3, nGmm = 0.1, nCrss=5){

######################Codes for test run
#  library("dplyr")
#  library("e1071")
#  library("stringr")
#  library("ggpubr")

#  Image = TRUE
#  nCst = 3
#  nGmm = 0.1
#  nCrss=5


########################################
  sdate <- Sys.Date()

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
    file_path1 <- choose.files(caption = "Select one or more files to set the directory (degenerate condition) ")

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
  file_path2 <- choose.files(caption = "Select one or more files to set the directory (intact condition) ")

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

  if(Image == TRUE){
    AiES::Axdistmap(Mainf[1], Sub_Back = TRUE, Binary = FALSE, All_Features = FALSE, Type = "tiff")
    AiES::Axdistmap(Mainf[2], Sub_Back = TRUE, Binary = FALSE, All_Features = FALSE, Type = "tiff")
  }
  #####################################

  txt_Group <- c("Degenerate","Intact")
  Ftr.tmp <- c("Group","m.eccentricity","s.radius.sd","h.sva.s2","h.idm.s1","h.sen.s1","m.majoraxis")

###############Create training date####

  if (exists('Data_sh') == TRUE) invisible({rm(Data_sh);gc();gc()})
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


##########################################
###Match the sample size of the two groups
  group_counts <- table(Data_sh$Group)
  min_group <- names(which.min(group_counts))
  filtered_min <- Data_sh[Data_sh$Group == min_group, ]

  min_counts <- group_counts[names(which.min(group_counts))]
  max_group <- names(which.max(group_counts))
  nrm_rows <- group_counts[names(which.max(group_counts))]-group_counts[names(which.min(group_counts))]
  major_group <- Data_sh$Group == max_group
  rm_rows <- sample(which(major_group), nrm_rows)

  #filtered_max <- sample(Data_sh[Data_sh$Group == max_group, ], min_counts)
  Data_sht <- Data_sh[!(major_group & row.names(Data_sh) %in% rm_rows),]
##########################################





##################delete unnecssary columns
  data_svm <- Data_sht[ which(colnames(Data_sht) %in% Ftr.tmp)]
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
  ######Selecting the Directory1#######
  # file full path
  file_path3 <- choose.files(default = paste0(sdate,"_AxClassifer.svm") ,caption = "Save a SVM model", multi = FALSE)

  if (length(file_path3)==0) {
    cat("Canceled \n")
    return(NULL)##Error countermeasure code
  }else{
    # extract directory info from the selected file
    dir_info3 <- dir_info(file_path3)
    colnames(dir_info3) <- c("Dir_Name", "Level")
    Mainf[3] <- str_c(dir_info3$Dir_Name[1:(nrow(dir_info3)-1)], collapse = "/")
    ##########################export data file
    setwd(Mainf[3])
    write.table(data_svm, paste0(sdate,"_Extracted_data_for_ML.txt"),
                sep="\t",row.names=FALSE, quote=F, col.names=TRUE, append=FALSE)
    #return(data_sht)
    ##########################################
    save(SVM_model, file= file_path3)
  }





}###end of function








