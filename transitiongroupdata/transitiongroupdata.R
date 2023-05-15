################
#author: YANG JUNJIE 
#MAY 2022
#data files were exported to the specified directory
#transitiongrou(file, num, m, filpath)
#num: number of transitions allowed in a RT group
#m: product ion value
#filepath: address for storing results
################
#execute in cmd line:
##C:\Program Files\R\R-4.1.0\bin>Rscript G:\zhengjie_project\MRM_pipline\transitiongroupdata\transitiongroup.R "G:\zhengjie_project\MRM_pipline\transitiongroupdata\OHgroup_example.csv" "50" "171.3" "G:\zhengjie_project\MRM_pipline\transitiongroupdata"
#Rscript transitiongroup_DnsClandMPEA_final.R website_transitionsinput_demo.csv 50#替换的这个 171.3 ./
#####################################################################################
#loading libraries in these functions
packsneed <- c('rcdk','rcdklibs','dplyr',"rio", 'hash')
lapply(packsneed, require, character.only = TRUE)

#intermid function for transition grouping function
build_intermediate_xlsx = function(source_data, 
                                   table, 
                                   export_num, 
                                   num, #number of transitions per minutes
                                   dir_name,  #differentiate the DnsCl and MPEA # seperated folder
                                   pi,   #specify product ion 1 '171.1' or 136.1
                                   pri2, #specify precursor ion 2: 2 or 3
                                   pi2, #specify product ion 2 '173.1' or 139.1
                                   ce,
                                   result_dir_name) {
  #loop over every output files
  for (i in 1:as.numeric(export_num)){
    #output data format
    result = data_frame()
    names(result) = c("CASRN",'Name', "SMILES", "MW","RT", "RTgroup", "MSgroup")
    #form each result for each ms group
    for (line in filter(table, n >= i)$MSgroup) {
      result = rbind(result, slice(filter(source_data, MSgroup == line), i))
    }
    #output the final results
    result = arrange(result, RTgroup)
    export(result, paste(dir_name, "/result", i, ".xlsx", sep = ""), which = "sheet1")
  }
  get_random_result(dir_name,
                    num, 
                    pi,   #specify product ion 1
                    pri2, #specify precursor ion 2
                    pi2, #specify product ion 2 
                    ce,
                    result_dir_name)
}

get_random_result = function(dir_name, 
                             num,
                             pi,   #specify product ion 1
                             pri2, #specify precursor ion 2 
                             pi2, #specify product ion 2
                             ce,
                             result_dir_name) {
  i = 1
  for (file in list.files(dir_name)) {
    # print(paste("---processing", i, "files---"))
    # if (file == "statistics.xlsx") {
    #   next
    # }
    data = import(paste(dir_name, "/", file, sep = ""))
    # RTgroup counting
    temp_table = arrange(data, RTgroup) %>% count(RTgroup) %>% 
      filter(n > num) %>% mutate(left_count = n - num)
    if (dim(temp_table)[1] == 0) {
      # if no further match for the MS group, then continue the searching
      format_export_data(result_dir_name, 
                         i, 
                         data, 
                         pi,   #specify product ion 1
                         pri2, #specify precursor ion 2
                         pi2, #specify product ion 2 
                         ce,
                         "sheet1")
      i = i + 1
      next
    }
    # store the row from the data out of the RT groups into a new datafile
    left_data = data_frame()
    names(left_data) = c("CASRN", 'Name', "SMILES", "MW","RT", "RTgroup", "MSgroup")
    for (line in temp_table$RTgroup) {
      # take tables for splitting, eliminate those not wanted
      left_num = filter(temp_table, RTgroup==line)$left_count
      RT_data = filter(data, RTgroup == line)
      left_data = rbind(left_data, RT_data[sample(nrow(RT_data), left_num, 
                                                  replace = F),])
    }
    data = setdiff(data, left_data)
    # order data
    data = arrange(data, RTgroup)
    left_data = arrange(left_data, RTgroup)
    # export data table
    format_export_data(result_dir_name,
                       i, 
                       data,
                       pi,   #specify product ion 1
                       pri2, #specify precursor ion 2
                       pi2, #specify product ion 2 
                       ce,
                       "sheet1")
    format_export_data(result_dir_name, 
                       i, 
                       left_data,
                       pi,   #specify product ion 1
                       pri2, #specify precursor ion 2
                       pi2, #specify product ion 2
                       ce,
                       "sheet2")
    i = i + 1
  }
  # delet temp files
  unlink(dir_name, recursive = TRUE)
}

#output formal reports
#format export for light & heavy
format_export_data = function (path, 
                               i, 
                               data,
                               pi,   #specify product ion 1
                               pri2, #specify precursor ion 2 '+2 or +3
                               pi2, #specify product ion 2
                               ce,
                               sheet_name) {
  data_length = dim(data)[1]
  #add condition for the light and heavy transitions
  #for molecuels with DnsCl
  pri2 = as.numeric(pri2)
  output1 = data.frame(
    "Compound Group" = rep("neg", data_length),
    "Compound Name" = data$Name, #
    "ISTD?" = rep("FALSE", data_length),
    "Precursor Ion" = data$MW + 1.0073, #
    "MS1 Res" = rep("Unit", data_length),
    "Product Ion" = rep(pi, data_length), #
    "MS2 Res" = rep("Unit", data_length),
    "Primary" = rep("TRUE", data_length),
    "Trigger" = rep("FALSE", data_length) ,
    "Ret Time (min)" = data$RT, #
    "Delta Ret Time" = rep("2", data_length),
    "Fragmentor" = rep("166", data_length) ,
    "Collision Energy" = rep(ce, data_length), #
    "Cell Accelerator Voltage" = rep("4", data_length) ,
    "Polarity" = rep("Positive", data_length),
    "Trigger Entrance Delay (cycles)" = rep("0", data_length),
    "Trigger Delay (cycles)" = rep("0", data_length),
    "Trigger Window" = rep("0", data_length),
    "IsLogicEnabled" = rep("FALSE", data_length),
    "Trigger Logic Flag" = rep("AND", data_length),
    "Trigger Ratio" = rep("1", data_length),
    "Trigger Ratio Window" = rep("1", data_length),
    "Ignore MRM" = rep("FALSE", data_length), check.names = F
  )
  #duplicate heavy CIL
  output2 = data.frame(
    "Compound Group" = rep("neg", data_length),
    "Compound Name" = data$Name, #
    "ISTD?" = rep("FALSE", data_length),
    "Precursor Ion" = data$MW + pri2 + 1.0073, #
    "MS1 Res" = rep("Unit", data_length),
    "Product Ion" = rep(pi2, data_length), #
    "MS2 Res" = rep("Unit", data_length),
    "Primary" = rep("TRUE", data_length),
    "Trigger" = rep("FALSE", data_length) ,
    "Ret Time (min)" = data$RT, #
    "Delta Ret Time" = rep("2", data_length),
    "Fragmentor" = rep("166", data_length) ,
    "Collision Energy" = rep(ce, data_length), #
    "Cell Accelerator Voltage" = rep("4", data_length) ,
    "Polarity" = rep("Positive", data_length),
    "Trigger Entrance Delay (cycles)" = rep("0", data_length),
    "Trigger Delay (cycles)" = rep("0", data_length),
    "Trigger Window" = rep("0", data_length),
    "IsLogicEnabled" = rep("FALSE", data_length),
    "Trigger Logic Flag" = rep("AND", data_length),
    "Trigger Ratio" = rep("1", data_length),
    "Trigger Ratio Window" = rep("1", data_length),
    "Ignore MRM" = rep("FALSE", data_length), check.names = F
  )
  #combine two output and export
  output = rbind(output1,output2)
  output = output %>% 
    arrange(output[,2])
  #output final files with experiment settings
  export(output, paste(path, "/result", i, ".xlsx", sep = ""), 
         which = sheet_name)
  #output data frame for step 4 peak analysis
  methodop1 = data.frame(
    "name" = data$Name,
    "RT" = data$RT,
    'MS1' = data$MW + 1.0073,
    'ID' = paste0(data$Name,"12"),
    'ID2' = data$Name
  )
  methodop2 = data.frame(
    "name" = data$Name,
    "RT" = data$RT,
    'MS1' = data$MW + pri2 + 1.0073,
    'ID' = paste0(data$Name,"13"),
    'ID2' = data$Name
  )
  step4method = rbind(methodop1, methodop2)
  step4method = step4method %>% 
    arrange(step4method[,1])
  export(step4method, paste(path,'/resultforstep4', i , '.xlsx', sep =''),
         which = sheet_name)
}

#parameters: num = number of transitions allowed permin, ce=collision energy, pi,pi2= production ion for DnsCl and MPEA, 
#pri2= mass increment for each precursor ion
transitiongroup = function(input, 
                           num, 
                           filpath,
                           pi,   #specify product ion 1
                           pri2, #specify precursor ion 2
                           pi2, #specify product ion 2 
                           ce) {
  # read data
  source_data = input
  # round number upwords, get rounded integer of MS1
  source_data = mutate(source_data, RTgroup=ceiling(RT), 
                       MSgroup=sprintf("%0.2f", ceiling(MW * 10) / 10))
  #count the number of compounds with same MS1 mass
  table = count(source_data, MSgroup)
  #set the number of transitions per minute into num/2, due to the heavy/light transitions, need to add heavy precursor and product in the row
  num = as.numeric(num)/2
  #grouping compounds into n files (n = count of compounds with same MS1 mass)
  export_num = slice(table, which.max(table$n))
  #find the folder for storing output files
  time_now = round(as.numeric(as.POSIXct(Sys.time())))
  dir_name = paste(filpath,"/temp",sep = "")
  dir.create(dir_name)
  #the folder for output result
  result_dir_name = paste(filpath,"/", "result_",sep = "")
  dir.create(result_dir_name)
  #output MSgroup and RT group data
  # export(table, paste(result_dir_name, "/statistics.xlsx", sep = ""))
  build_intermediate_xlsx(source_data, 
                          table, 
                          export_num$n, 
                          num, 
                          dir_name, 
                          pi,   #specify product ion 1
                          pri2, #specify precursor ion 2
                          pi2, #specify product ion 2 
                          ce,
                          result_dir_name)
}

############################################final functions####################################
#execute function for transitions grouping
transitionsplit = function(
                      inputfile, 
                      num, 
                      filpath,
                      output){
  datinput = read.csv(inputfile)

  
  #############function for splitting MPEA and DnsCl
  #sort the data by DnsCl and MPEA
   #updated on Oct 2022, need to update the query fragment
  query1 =  c('c1(ccc(cc1)-[#7](-[#6])-[#6])') ######DnsCl for OH model
  query2 = c('c1ccc(cc1)-[#6]-[#6]-[#7](-[#6])-[#6](-[#6])=[#8]') #####MPEA for COOH model
  
  #####Oct2022 if dnscl_smiles or mpea_smiles is empty, there will be error
  if (sum(!is.na(datinput['DnsCl_SMILES']))> 0){
    datinput2 = datinput %>% 
      rowwise() %>%
      mutate(mols1 = parse.smiles(DnsCl_SMILES))
    datinput2 = datinput2 %>% mutate(DnsCl = case_when(
      rcdk::matches(query1,mols1) ~ 'y',
      TRUE ~ 'n'))
  } else{
    datinput2 = datinput %>% mutate(DnsCl = 'n')
  }
  
  if (sum(!is.na(datinput['MPEA_SMILES'])) > 0){
    datinput2 = datinput2 %>% 
      rowwise() %>%
      mutate(mols2 = parse.smiles(MPEA_SMILES)) 
    datinput2 = datinput2 %>% mutate(MPEA= case_when(
      rcdk::matches(query2,mols2) ~ 'y',
      TRUE ~ 'n'))
  }else{
    datinput2 = datinput2 %>% mutate(MPEA ='n')
  }
  #separate rows with DnsCl and MPEA
  datinput_1 = datinput2 %>% filter(DnsCl =='y')
  datinput_2 = datinput2 %>% filter(MPEA == 'y')
  #select columns and renames column name, set the same names as the transitiongroup functions
  datinput_1 = datinput_1 %>% select(Name,DnsCl_SMILES,DnsCl_mass, DnsCl_RT)
  datinput_2 = datinput_2 %>% select(Name,MPEA_SMILES,MPEA_mass, MPEA_RT)
  colnames(datinput_1) = c('Name','SMILES','MW','RT')
  colnames(datinput_2) = c('Name','SMILES','MW','RT')
  datinput_1 = as.data.frame(datinput_1)  #convert tibble to dataframes for post-processing
  datinput_2 = as.data.frame(datinput_2)

  dir_name = paste(filpath,"/",output, sep = "")
  dir.create(dir_name)
  #split MRM transitions for DnsCl and MPEA
   if (nrow(datinput_1)> 1){
    #path1 = paste0(dir_name,'/DnsCl', format(Sys.time(), "%F %H-%M"))
    path1 = paste0(dir_name, '/DnsCl')
    transitiongroup(datinput_1,
                    num,
                    path1,
                    pi = 171.1,   #specify product ion 1
                    pri2 = 2, #specify precursor ion 2
                    pi2 = 173.1, #specify product ion 2
                    ce = 35)
  }
  
  if (nrow(datinput_2)> 1){
    #path2 = paste0(dir_name, '/MPEA', format(Sys.time(), "%F %H-%M"))
    path2 = paste0(dir_name, '/MPEA')
    dir.create(path2)
    transitiongroup(datinput_2,
                    num,
                    path2,
                    pi = 136.1,   #specify product ion 1
                    pri2 = 3, #specify precursor ion 2
                    pi2 = 139.1, #specify product ion 2
                    ce = 40)
  }
}
 

###################################
#take command line string as inputs
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
num = args[2]
filpath = args[3]
output = args[4]

#execute function
transitionsplit(file,num,filpath,output)

#create data table for direct input the step 4


#####################final grouping results
####function test
# setwd('D:/zhengjie_project/MRM_pipline/transitiongroupdata')
# demo = 'website_transitionsinput_demo.csv'
# num = 50
# filpath ='./'

args <- commandArgs(trailingOnly = TRUE)
file = args[1]
num = args[2]
filpath = args[3]
output = args[4]
transitionsplit(file,
                num,
                filpath,
                output)

