#author: YANG JUNJIE date: May 2022
#datafile generated in the filepath
#how to run in cmd line
#C:\Program Files\R\R-4.1.0\bin>Rscript G:\zhengjie_project\MRM_pipline\predictionRT_script\predictionRT.R "G:\zhengjie_project\MRM_pipline\predictionRT_script\OH_demoRT.csv" "OH" "G:/zhengjie_project/MRM_pipline" "result"

# Rscript predictionRT_DnsClandMPEA_debug.R OH_demoRT_test.csv OH ./ result123
#######################################################################################
#------------------------rt prediction module------import rt model--------------------
#loading packages
#load all packages at once
# setwd('D:/zhengjie_project/MRM_pipline/predictionRT')
# getwd()


##########retention time prediction regardless of light and heavy CIL
#newfunction for RT predictionp
predictRT = function(input,
                     filpath){
  datinput = input
  ####splitting MPEA and DnsCl
  #sort the data by DnsCl and MPEA
  query1 =  c('CN(C)c(ccc1)c2c1c(S(=O)(OC)=O)ccc2') #DnsCl for OH model
  query2 = c('CN(C)CCc1ccccc1') #MPEA for COOH model
  datinput = datinput %>% 
    rowwise() %>%
    mutate(mols1 = parse.smiles(DnsCl_SMILES))%>%
    rowwise() %>%
    mutate(mols2 = parse.smiles(MPEA_SMILES))
  
  datinput2 = datinput %>% mutate(DnsCl= case_when(
    rcdk::matches(query1,mols1) ~ 'y',
    TRUE ~ 'n'
  )) %>% mutate(MPEA = case_when(
    rcdk::matches(query2,mols2) ~ 'y',
    TRUE ~ 'n'
  ))
  #separate rows with DnsCl and MPEA
  datinput_1 = datinput2 %>% 
    filter(DnsCl =='y') %>%
    select(-c(SMILES))      #remove original smiles
  datinput_1 = datinput_1 %>%
    rename(SMILES = DnsCl_SMILES)  #rename derivatized smiles as SMILES for RTprediction
  datinput_2 = datinput2 %>% 
    filter(MPEA == 'y') %>%
    select(-c(SMILES)) 
  datinput_2 = datinput_2 %>%
    rename(SMILES = MPEA_SMILES)

  #predict and combine all dataframe and output
  predrt1 = as.data.frame(predictionRT(datinput_1,
                         'DnsCl', 
                         filpath))
  predrt2 = as.data.frame(predictionRT(datinput_2, 
                         'MPEA', 
                         filpath))
  print(dim(predrt1))
  print(dim(predrt2))
  #coombine all columns with the same smiles and add another column
  datinput = datinput %>% 
    mutate_at(c('DnsCl_SMILES','MPEA_SMILES'), ~na_if(., '')) %>%
    select(-contains('mol')) #drop mols1 and mols2
  datinput = datinput %>%
    left_join(predrt1, by=c('DnsCl_SMILES'='DnsCl_SMILES'))
  print(dim(datinput))
  newoutput = datinput %>%
    left_join(predrt2, by=c('MPEA_SMILES'='MPEA_SMILES')) %>%
    distinct()
  print(dim(newoutput))
  return(newoutput)  #final output
}


##a wrapper function for prediction RT to be launched by CML
predictionRT <- function(input,b,filpath){
  #RUN FUNC
  newsmilst <- input$SMILES
  getdesc_nopp <- function(input2){
    mols <- parse.smiles(input2)
    descNames <- unique(unlist(sapply(get.desc.categories(),get.desc.names)))
    ## todo
    print(descNames)
    descNames <- descNames[descNames!= 'org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor']
    descs_tot <- data.frame()
    descs_temp <- eval.desc(mols,descNames,verbose = FALSE)
    return(descs_temp)
  }
  descs <- getdesc_nopp(newsmilst) #get mds without any data imputation
  #find a way to load colnames data, settings,and model
  #define functions for molecular descriptors without preprocess
  load.Rdata(paste0(filpath,'/',b,'_descs.RData'),'descname')
  load.Rdata(paste0(filpath,'/',b,'_normsettings.RData'),'settings')
  load.Rdata(paste0(filpath,'/',b,'_model.RData'),'predmodel')
  descsdata <- descs[,descname] #select mds to prepare same dimensions for normalzation settings
  descsdata <- predict(settings,descsdata) #get normalized mds for modeling
  predictionrt <- predict(predmodel, newdata = descsdata) #model select the mds as it trained
  output = data.frame(input$SMILES, predictionrt)
  # output <- cbind(newsmilst, predictionrt)
  smicolname = paste0(b,'_SMILES')
  RTcolname = paste0(b,"_RT")
  colnames(output) <- c(smicolname,RTcolname)
  print(head(output))
  # output %>% export(paste0(filpath,"/",result,".csv")) #csv import
  return(output)
}



packsneed <- c('rcdk','rcdklibs','randomForest','leaps','caret','corrplot',
               'mlr','dplyr','Metrics','ggpubr','ggplot2','miceadds','rio')
lapply(packsneed, require, character.only = TRUE)

#take command line string as inputs
args <- commandArgs(trailingOnly = TRUE)
input = args[1]
b=args[2]
filpath = args[3]
result = args[4]

newinput <- read.csv(input)
# predictionRT(newinput, b, filpath,result)
#filepath change to ./

#demo data
#demoRT = read.csv('website_input_for_RT.csv')
setwd("/home/ubuntu/cooh_server/predictionRT")
demoRT_result = predictRT(newinput, getwd())
setwd("/home/ubuntu/cooh_server/static")
demoRT_result %>% export(paste0(result,".csv"))
print(demoRT_result)
