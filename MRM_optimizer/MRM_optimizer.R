#author: YANGJUNJIE 
#date:Sep 2022
#rtinput:
#mspath:
#respath:
#execute in cmd line:
#C:\Program Files\R\R-4.1.0\bin>Rscript G:\zhengjie_project\MRM_pipline\peakanalysis\peakanalysis.R "G:\zhengjie_project\MRM_pipline\peakanalysis\method2_MS1RT.csv" "G:\zhengjie_project\MRM_pipline\peakanalysis" 0.01 30 0.1 1.5 300
#cmd line in linux:
#C:\Program Files\R\R-4.1.0\bin>Rscript .\peakanalysis.R .\method2_MS1RT.csv  .\peakanalysis 0.01 30 0.1 1.5 300

#Rscript ./peakanalysis_final.R method2_MS1RT.csv ./ 0.01 30 0.1 1.5##替换这个 300
##############test on workstation folder
# rtinput = 'method2_MS1RT.csv'
# mspath = 'G:/zhengjie_project/MRM_pipline/peakanalysis'
# pwid1 = 0.01
# pwid2 = 30
# snthr = 0.1
# delrt = 1.5
# intenmin = 300
##############test on PC folder
# pwid1 = 0.01
# pwid2 = 30
# snthr = 0.1
# delrt = 1.5
# intenmin = 300
# specpath = './peakanalysis/rawMSdata2/'
# rtinput = "./peakanalysis/method2_MS1RT.csv"
# peakanalysis(rtinput=rtinput, 
#              mspath=specpath,
#              pwid1=pwid1,
#              pwid2=pwid2,
#              snthr=snthr,
#              delrt=delrt,
#              intenmin=intenmin)

# #####################example of data input for this function#################
# args <- commandArgs(trailingOnly = TRUE)
# rtinput = args[1]
# mspath = args[2]
# pwid1 = as.numeric(args[3])  #how to record this parameters #direct input c(0.01,30)?
# pwid2 = as.numeric(args[4])
# snthr = as.numeric(args[5])
# delrt = as.numeric(args[6])
# intenmin = as.numeric(args[7])
# # resname = args[8]
# peakanalysis(rtinput,mspath,pwid1,pwid2,snthr,delrt,intenmin)

#############################################FUNCTION#########################################################
#loading packages
packsneed <- c('xcms','magrittr','MSnbase','dplyr','tidyr','ggplot2','tidyverse','ggpubr',"ggrepel","rio", 'ggfortify')
lapply(packsneed, require, character.only = TRUE)

#function for screening peaks with both light and heavy signals
#cleaning peaks with both heavy and light signals
data_process = function (raw_data, wrong_data, output) {
  # 维护一个遍历数据
  category_list = unique(raw_data$ID2)
  for (category in category_list) {
    # 查出同一组的行
    data = filter(raw_data, ID2==category)
    # 首先判断是否合法数据 12 13个数
    data_12_line_num = c()
    data_13_line_num = c()
    for (i in 1:nrow(data)) {
      row_data = data[i, ]
      if (grepl(12, row_data$ID)) {
        data_12_line_num = c(data_12_line_num, i)
      }
      if (grepl(13, row_data$ID)) {
        data_13_line_num = c(data_13_line_num, i)
      }
    }
    count_12 = length(data_12_line_num)
    count_13 = length(data_13_line_num)
    if (count_12 == 0 || count_13 == 0) {
      # 有任何一个为0 视为错误数据
      wrong_data = rbind(wrong_data, data)
    } else if (count_12 > count_13) { # 正确数据中选出小的那个作为标准
      ret = data_filter(data, data_12_line_num, data_13_line_num, 
                               wrong_data, output)
      wrong_data = ret$w
      output = ret$o
    } else {
      ret = data_filter(data, data_13_line_num, data_12_line_num, 
                               wrong_data, output)
      wrong_data = ret$w
      output = ret$o
    }
  }
  return(list(w=wrong_data, o=output))
}

#function for tag the light and heavy isotope amount the data, add the 4th ID for further matching
#data is the peak table
data_filter = function (data, num_1, num_2, wrong, output) {
  # 以num_2为基准
  num = 1
  for(i in num_2) { # 分组
    chosen_row = data[i, ]
    quantity = 0
    for (j in num_1) {
      compare_row = data[j, ]
      # 判断条件
      deltart2 =abs(chosen_row$rt-compare_row$rt)
      deltamaxo = chosen_row$maxo/compare_row$maxo
      if (deltart2 <= 0.1 && deltamaxo >= 0.25 && deltamaxo <= 4) {
        # 同时满足这三条留下
        # 暂时不考虑一个对多个符合的情况
        quantity = quantity + 1
        #add a new 4th name for later cross matching
        compare_row$ID3 = paste(compare_row$ID2, "-", num, sep="")
        chosen_row$ID3 = paste(compare_row$ID2, "-", num, sep="")
        output = rbind(output, compare_row)
        output = rbind(output, chosen_row)
      }
    }
    num = num + 1
    if (quantity != 1) {
      wrong = rbind(wrong, data)
      return(list(w=wrong, o=output))
    }
  }
  return(list(w=wrong, o=output))
} 

#get_category: get clean heavy and light peak list for each sample group
#input: raw peak table and the number of sample groups
#num: number of groups, since column value is the name of sample group
#output: clean peak table with ID3
get_category = function (data, snum) {
  #define wrong data
  wrong_data = data.frame()
  #define right data
  output = data.frame()
  for(num in 1:snum) {
    w_data = data.frame()
    o_data = data.frame()
    datalist = filter(data, column == num)
    #data_process: 
    ret = data_process(datalist, w_data, o_data)
    wrong_data = rbind(wrong_data, ret$w)
    output = rbind(output, ret$o)
  }
  return(list(w=wrong_data, o=output))
}

#define the peak detection function and generate the raw peak table
peaksdetect <- function(input,rtdata,pwid1,pwid2,snthr, delrt, 
                        intenmin,
                        sample){
  dmrm <- readSRMData(input)
  cwp <- CentWaveParam(peakwidth = c(pwid1,pwid2), snthresh = snthr)#snthresh ratio adjust to 10, peakwidth to 0.01-30
  grouppeaks <- findChromPeaks(dmrm,cwp)
  mspeaks <- chromPeaks(grouppeaks)
  rawpeaks <- as.data.frame(mspeaks)
  rawpeaks['Mz1'] <- precursorMz(grouppeaks)[mspeaks[, "row"], ][,1]
  rawpeaks['Mz2'] <- productMz(grouppeaks)[mspeaks[, "row"], ][,1]
  peaktable <- rawpeaks %>%
    left_join(rtdata, by='row')
  #peak screening
  #using the Predicted RT for the filtering  #rt <=1.5
  #using the intensity threshold filtering   maxo >= intenmin
  peaktable <- peaktable %>%
    mutate('deltart'= abs(rt-predrt)) %>%
    filter(deltart <= delrt & maxo >= intenmin) %>%
    left_join(sample, by = c('column'='rowID')) %>%   ##combine sample information with peak table
    select(rt,maxo,row,column,Mz1,Mz2,ID, deltart,ID2, sample_name, sample_group)
  return(peaktable)
}

####FUNC for refine peak table for statistics analysis 
#input: peak table & the names of sample groups
#input peak table with ID3
export_func = function (data,sname) {
  finaldat = data.frame()
  # peak = peaklist #unique peak id
  peak = unique(data$ID3)
  snum = length(sname)
  #add sname to sample_name
  sample_name = c("ID3",sname)
  #get total num of sample group
  totnum = length(sname)
  #collect peaks according to peakid
  for (id3 in peak) {
    temprow = c(id3)
    #find out the bad rows in each treatment group
    for (i in 1:totnum) {
      if (nrow(filter(data, column == i & ID3 == id3)) == 0) {
        temprow[i+1] = NA
        next
      }
      temprow[i+1] = filter(data, column == i & ID3 == id3)$maxo
    }
    finaldat = rbind(finaldat, temprow)
  }
  names(finaldat) = sample_name
  return(finaldat)
  #what is finaldat format
}

peaksep = function (data,snum) {
  final = data.frame()
  peak = unique(data$ID3)
  peak2 = c()
  for (id3 in peak) {
    nacount1 = 0
    # nacount2 = 0
    #find out the bad rows in all groups (each treatment group), how to remove in each group
    for (i in 1:snum) {
      if (nrow(filter(data, column == i & ID3 == id3)) == 0) {
        nacount1 = nacount1 + 1 
        next
      }
    }
    if (nacount1 <= 1){
      peak2 = c(peak2,id3)
    } else{
      next
    }
  }
  return(peak2)
}

#refine data format to take multiple samples
peakrefine <- function(input,sname){
  group_num = length(sname)
  res = get_category(input, group_num)
  output = res$o
  rownames(output) = seq(1, nrow(output), 1)
  output = filter(output, grepl("12", ID)) # only consider isotope-12 peaks
  retain_peak = peaksep(output, group_num)
  #further remove peak id that can not be detected in all groups
  final_output = filter(output,ID3%in%retain_peak)
  # export_func: input: filtered peak table with ID3 & the name of sample group
  final_table = export_func(final_output,sname) #group1 and group2 as parameters
  rownames(final_table) = final_table$ID3
  #analyzetable: final output peak table
  analyzetable = subset(final_table,select = -c(ID3))
  return(analyzetable)
}

#ttest analysis   #maybe take one-way anova
ttestpeaks <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y,
                   alternative=c("two.sided"),
                   paired=FALSE,var.equal=FALSE)
  results$p.value
  return(results$p.value)
}

#define the function for differential express and graphs output
diffexp = function(input,grp1, grp2,output){
  pvalue <- apply(input, 1, ttestpeaks, grp1 = grp1, grp2 = grp2)
  #replace nas in each group with group means
  Log2FC = apply(input[,grp2],1,mean) - apply(input[,grp1],1,mean)
  res1 <- cbind(Log2FC,pvalue)
  res1 <- as.data.frame(res1)
  res1$diffexpressed <- "NO" #add a column of NAs #add label for down and up regulation
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  res1$diffexpressed[res1$Log2FC > 0.6 & res1$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res1$diffexpressed[res1$Log2FC < -0.6 & res1$pvalue < 0.05] <- "DOWN"
  res1$delabel <- NA
  res1$symbol <- rownames(res1)
  res1$delabel[res1$diffexpressed != "NO"] <- res1$symbol[res1$diffexpressed != "NO"]
  #volcanoplot
p <- ggplot(data=res1, aes(x=Log2FC, y=-log10(pvalue),fill= diffexpressed, label=delabel)) +
    geom_point(size=3, shape=21) + 
    scale_color_manual(values=c("#999999", "#0066CC", "#56B4E9"))+
    theme(text = element_text(size=15, family = 'Arial',face='bold'),
          axis.text.x = element_text(size=15, family = 'Arial',face='bold',colour = "black"),
          axis.text.y = element_text(size=15, family = 'Arial',face='bold',colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1), # Black border around the plot area.
          axis.ticks = element_line(colour = "black", size = 1.5), # Style of x and y ticks.
    ) +
    geom_hline(yintercept = 1.30103, colour = "black", linetype = "dashed", size = 0.75) + # Horizontal significance cut-off line.
    geom_vline(xintercept = 0.584963, colour = "black", linetype = "dashed", size = 0.75) + 
    annotate("rect", xmin = 0.584963, xmax = Inf, ymin = 1.30103, ymax = Inf, alpha = .2) # Shade plot subregion.
  # print(p) #return the plot
  name = paste0(output,".png")
  setwd("/home/ubuntu/cooh_server/static")
  ggsave(filename = name, plot = p, bg='white', units="in", width = 8, height=6, dpi=400)
  setwd("/home/ubuntu/cooh_server")
  return(res1)
}

######################################################
#function for peaks analysis #peak and data extraction
peakanalysis <- function(rtinput, 
                         mspath,
                         pwid1,
                         pwid2,
                         snthr,
                         delrt,
                         intenmin,
                         output)
  {
  MSRT = read.csv(paste0("./",rtinput))
  rt = data.frame("row"=1:length(MSRT$RT),"ms1"=MSRT$MS1, "predrt"=MSRT$RT, "ID"=MSRT$ID,"ID2" = MSRT$ID2)
  sample = dir(path = mspath, pattern = 'samplename.csv$', recursive = TRUE, full.names = TRUE)
  print(sample)
  sampleinfo = read.csv(paste0(mspath,"/samplename.csv"))
  sample_name = sampleinfo[,'sample_name']
  sample_group = sampleinfo[,'sample_group'] 
  print(mspath)
  print(length(unique(sample_group)))
  filedat = list.files(path = mspath, pattern = ".mzML$", full.names = TRUE)
  print(filedat)
  #generate peaktable using dplyr pipline
  peaktable = peaksdetect(filedat, 
                          rt, 
                          pwid1,
                          pwid2,
                          snthr,
                          delrt,
                          intenmin,
                          sampleinfo) #these function need to be refined for their result address
  #clean up peak table with more filters
  refinetable = peakrefine(peaktable,
                            sample_name)
  #perform peak analysis
  #select the analysis methods
  if (length(unique(sample_group))>2){
    print('using ANOVA for analysis')
    #convert all list to numeric number 
    refine_peak = as.data.frame(apply(refinetable, 2, as.numeric))
    #add rownames as ID
    refine_peak = cbind(ID = rownames(refinetable), refine_peak)
    #change row name to sample group id, pivot the table in order to making anova analysis
    peak_pivot = refine_peak %>%
      pivot_longer(cols = 2:ncol(.)) %>%
      pivot_wider(names_from = ID, values_from = value) %>%
      separate(name, into = c("trt"), sep = "_")  #rename the rowname for easier analysis
    #take log2 value for all 
    peak_pivot2 = cbind(peak_pivot[,1],apply(peak_pivot[,2:ncol(peak_pivot)], 2, log2))
    #replace NA with 0 signal 
    peak_pivot2[is.na(peak_pivot2)] = 0
    #perform ANOVA analysis 
    fit_aov <- function(col) {
      anovas = aov(col ~ trt, data = peak_pivot2)
      result =TukeyHSD(anovas)
      resultdf = data.frame(result$'trt')
      resultdf = resultdf[,c('diff','p.adj')]
      return(resultdf)
    }
    #generate list-like results
    anovas1 <- map(peak_pivot2[,2:ncol(peak_pivot2)], fit_aov)
    #output results in a csv
    anovas1df = as.data.frame(do.call(rbind, anovas1))
    write.csv(anovas1df, file='./anovas_result.csv')
  } else if (length(unique(sample_group) <= 2)){
    print('using Student T test for analysis')
    #perform Student T test analysis
    #take log value
    #replace NA with 0 signal 
    #ttest analysis   #maybe take one-way anova
    testtable = as.data.frame(apply(refinetable, 2, as.numeric))
    row.names(testtable) = row.names(refinetable)
    testtable = log2(testtable)
    testtable[is.na(testtable)]=0
    #get grp1 and grp2
    grp1 = colnames(testtable)[grep('ctr',colnames(testtable))]
    grp2 = colnames(testtable)[grep('tre',colnames(testtable))]
    diffres = diffexp(testtable, grp1,grp2,output)
    diffres = diffres %>%
      select(-delabel)
    setwd("/home/ubuntu/cooh_server/static")
    diffres %>% export(paste0(output,".csv"))
    setwd("/home/ubuntu/cooh_server/static")
  } else{
    print('data is not complete')
  }
}

##################################################code test############################################
# pwid1 = 0.01
# pwid2 = 30
# snthr = 0.1
# delrt = 1.5
# intenmin = 300
# specpath = './peakanalysis/rawMSdata2/'
# rtinput = "./peakanalysis/method2_MS1RT.csv"
# peakanalysis(rtinput=rtinput, 
#              mspath=specpath,
#              pwid1=pwid1,
#              pwid2=pwid2,
#              snthr=snthr,
#              delrt=delrt,
#              intenmin=intenmin)

#####################example of data input for this function#################
args <- commandArgs(trailingOnly = TRUE)
rtinput = args[1]
mspath = args[2]
pwid1 = as.numeric(args[3])  #how to record this parameters #direct input c(0.01,30)?
pwid2 = as.numeric(args[4])
snthr = as.numeric(args[5])
delrt = as.numeric(args[6])
intenmin = as.numeric(args[7])
output = args[8]
# resname = args[8]
peakanalysis(rtinput,mspath,pwid1,pwid2,snthr,delrt,intenmin,output)