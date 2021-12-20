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
model.sel.wt.ko = function(data.wt.ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18), period=24)
{
  #index =138; data.wt = as.numeric(aa[index, c(1:16)]);data.com = as.numeric(aa[index, grep('Cry.KO', colnames(nuclear))]);data.wt.ko =  c(data.wt, data.com); time.wt=c(0:15)*3; time.ko=c(0,6,12,18);period=24;
  wt = as.numeric(data.wt.ko[1:16])
  corr.ko = correlation.wt.ko(data.wt.ko)
  kk = which(!is.na(wt)==TRUE)
  wt = wt[kk]
  time.wt = time.wt[kk]
  ko = as.numeric(data.wt.ko[17:20])
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
nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')
#load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Table_Nuclear_mRNA_Total_Nascent.Rdata')
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
load(file='Tables_DATA/Table_Nuclear_total_mRNA_pol2.Rdata')
table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v3.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
source('functions_nuclear.R')
#kk = which(table.sx$Gene.names=='Sirt7')
res.ko = c()
for(n in 1:nrow(table.sx))
{
  index = table.sx[n, 1]
  data.wt = as.numeric(nuclear[index, c(1:16)]);
  res.ko = rbind(res.ko, model.sel.wt.ko(c(data.wt, as.numeric(nuclear[index, grep('Cry.KO', colnames(nuclear))]))))
  #cor.wtko = c(cor.wtko, correlation.wt.ko(data.wt, data.com))
  #diff.mean = c(diff.mean, (mean(data.com[which(!is.na(data.com)==TRUE)])-mean(data.wt[which(!is.na(data.wt)==TRUE)])))
  #correlation = c(correlation, correlation.wt.ko(c(data.wt, data.com)))
}

table.sx$prob.rhythmic.same.parameters.wt.ko = res.ko[,1]





