##########################################################################
##########################################################################
# Project:
# Script purpose: comapre the rhythmicity bewteen WT and KO for Daniel  
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Dec 20 12:02:39 2021
##########################################################################
##########################################################################

##########################################
# the main function used for model comparison 
##########################################
correlation.wt.ko = function(wt, ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18), period = 24)
{
  cat('calculate the correcitons between wt and ko for shared time points \n')
  
  # search for overlapping time points between wt and ko, if found, average of replicates were calculated
  wt.mean = c()
  tt.wt = time.wt%%period
  for(t in time.ko)
  {
    wt.mean = c(wt.mean, mean(wt[which(tt.wt == t)]))  
  }
  
  # wt1 = data.wt.ko[c(1,3,5,7)]
  # wt2 = data.wt.ko[c(9,11,13,15)]
  # ko = data.wt.ko[17:20]
  # wt = c()
  # for(n in 1:4)
  # {
  #   wt12 = c(wt1[n], wt2[n])
  #   wt = c(wt, mean(wt12[which(!is.na(wt12))]))
  # }
  # 
  
  if(all(!is.na(wt.mean)) && all(!is.na(ko)))
  {
    cat('correlation is cacluated only if all wt and ko at shared time points are not NA \n')
    return(cor(wt.mean, ko))
  }else{
    return(NA)
  }
}

model.sel.wt.ko = function(data_wt, data_ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18), period=24)
{
  #index =138; data.wt = as.numeric(aa[index, c(1:16)]);data.com = as.numeric(aa[index, grep('Cry.KO', colnames(nuclear))]);data.wt.ko =  c(data.wt, data.com); time.wt=c(0:15)*3; time.ko=c(0,6,12,18);period=24;
  wt = as.numeric(data_wt)
  ko = as.numeric(data_ko)
  corr.ko = correlation.wt.ko(wt, ko, time.wt, time.ko)
  
  kk = which(!is.na(wt)==TRUE)
  wt = wt[kk]
  time.wt = time.wt[kk]
  
  jj = which(!is.na(ko)==TRUE)
  ko = ko[jj]
  time.ko = time.ko[jj]
  nb.ko = length(ko)
  
  if(length(wt)<=8|length(ko)<=3)
  {
    prob.wt.ko = c(NA, NA, NA)
  }else{
    ### fitting wt.ko pool with rhythmic parameters
    wt.ko = c(wt, ko)
    time.wt.ko = c(time.wt, time.ko)
    c=cos(2*pi*time.wt.ko/period)
    s=sin(2*pi*time.wt.ko/period)
    fit = lm(wt.ko~c+s)
    rss = sum(fit$residuals^2)
    
    ###fitting wt with rhythmic parameter
    c1=cos(2*pi*time.wt/period)
    s1=sin(2*pi*time.wt/period)
    fit1 = lm(wt~c1+s1)
    rss1 = sum(fit1$residuals^2)
    
    ### fitting ko with rhythmic parameters
    c2=cos(2*pi*time.ko/period)
    s2=sin(2*pi*time.ko/period)
    fit2 = lm(ko~c2+s2)
    rss2 = sum(fit2$residuals^2)
    
    ### fitting ko with flat parameters
    rss3 = sum((ko-mean(ko))^2)
    
    ## model 1: rhythmic with same parameters
    n = length(wt.ko)
    BIC1 = n*log(rss/n) + 3*log(n);
    BIC2 = n*log((rss1+rss2)/n) + 6*log(n);
    BIC3 = n*log((rss1+rss3)/n) + 4*log(n);
    
    BIC = c(BIC1, BIC2, BIC3)
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    
    prob.wt.ko = prob.model
  }
  
  #print(pval.wt.ko)
  names(prob.wt.ko) = c('prob.M1', 'prob.M2', 'prob.M3')
  return(c(prob.wt.ko, nb.ko=nb.ko, corr.ko =corr.ko))
}

##########################################
# example how to use the main function 
##########################################
dataDir = '../data/nuclearProt/Tables_DATA/' # specific the folder of example table

nuclear = read.table(paste0(dataDir, 'nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt'), 
                     sep='\t', header=TRUE, as.is=c(17:20))

res = c()
for(n in 1:nrow(nuclear))
{
  # n = 1
  data.wt = nuclear[n, c(1:16)]
  data.ko = nuclear[n, grep('Cry.KO', colnames(nuclear))]
  
  res = rbind(res, model.sel.wt.ko(data_wt = data.wt, # wt data
                                   data_ko = data.ko, # ko data
                                   time.wt=c(0:15)*3, # time point of wt
                                   time.ko=c(0,6,12,18), # ko time points
                                   period=24) # period, 24h normally
              ) 
  
}

