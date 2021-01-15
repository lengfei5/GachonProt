#################################################################################### function for protein complex analysis
library(Hmisc)
library(CircStats)
cutoff.rhy = 0.01
cutoff.sim = 0.80

########################################################
########################################################
# Section : general utility functions
#  
########################################################
########################################################
first.upper = function(character)
{
    test = character
    test = tolower(test)
    test = c(toupper(unlist(strsplit(as.character(test),""))[1]), unlist(strsplit(as.character(test),""))[-1])
    test = paste(test, collapse = '')
    return(test)
}

## fisher method to combine p values
Fisher.test <- function(p) 
{
	Xsq <- -2*sum(log(p))
	p.val <- 1-pchisq(Xsq, df = 2*length(p))
	return(p.value = p.val)
}

#p <- c(.001, .2, .3)
#Fisher.test(p = p)

index.outliers = function(data.xx)
{
	c = 1.5
	Q1 = quantile(data.xx, 0.25,type=5)
	Q3 = quantile(data.xx, 0.75, type=5)
	IQD = Q3 - Q1
	lower = Q1 - c*IQD
	upper = Q3 + c*IQD
	index = which(data.xx<lower|data.xx>upper)
}


mutual.correlation = function(v1, v2, cutoff.nb.timepoints=6)
{
	t1 = which(!is.na(v1))
	t2 = which(!is.na(v2))
	t = intersect(t1, t2)
#return(abs(cor(v1[t], v2[t])))
	return(cor(v1[t], v2[t]))
#if(length(t)>=cutoff.nb.timepoints) {
#	return(corr(v1[t], v2[t]))
#}else{
#	return(NA)
#}
	
}
phase.difference = function(kk)
{
	phase.diff = c()
	ii = c(1:(length(kk)-1))
	for(i in ii)
	{
		jj = c(1:length(kk))
		jj = jj[which(jj>i)]
		for(j in jj)
		{
			diff = nuclear$phase[kk[j]]-nuclear$phase[kk[i]];
			if(diff>12) diff = diff-24;
			if(diff<(-12)) diff = diff +24;
			
			phase.diff = c(phase.diff, diff)
		}
	}
	return(phase.diff)
}

########################################################
########################################################
# Section : functions of mutatn vs wt 
# 
########################################################
########################################################
#### correlation between WT (16 time points) and KO (4 time points)
correlation.wt.ko = function(data.wt.ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18))
{
  wt1 = data.wt.ko[c(1,3,5,7)]
  wt2 = data.wt.ko[c(9,11,13,15)]
  ko = data.wt.ko[17:20]
  wt = c()
  for(n in 1:4)
  {
    wt12 = c(wt1[n], wt2[n])
    wt = c(wt, mean(wt12[which(!is.na(wt12))]))
  }
  if(all(!is.na(wt)) && all(!is.na(ko)))
  {
    return(cor(wt, ko))
  }else{
    return(NA)
  }
}
chow.test = function(data.wt.ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18), period=24)
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
  if(length(wt)<4|length(ko)<3)
  {
    pval.wt.ko = NA
  }else{
    ### fitting wt.ko pool
    wt.ko = c(wt, ko)
    time.wt.ko = c(time.wt, time.ko)
    c=cos(2*pi*time.wt.ko/period)
    s=sin(2*pi*time.wt.ko/period)
    fit = lm(wt.ko~c+s)
    ee = sum(fit$residuals^2)
    
    ###fitting wt
    c1=cos(2*pi*time.wt/period)
    s1=sin(2*pi*time.wt/period)
    fit1 = lm(wt~c1+s1)
    ee1 = sum(fit1$residuals^2)
    
    if(length(ko)>3)
    {
      ### fitting ko
      c2=cos(2*pi*time.ko/period)
      s2=sin(2*pi*time.ko/period)
      fit2 = lm(ko~c2+s2)
      ee2 = sum(fit2$residuals^2)
      
      F2 = (ee-ee1-ee2)/3/((ee1+ee2)/(length(wt)+length(ko)-2*3))
      pval.wt.ko = pf(F2, 3, (length(wt)+length(ko)-2*3), lower.tail = FALSE, log.p = FALSE)
    }else{
      F1 = (ee-ee1)/length(ko)/(ee1/(length(wt)-3))
      pval.wt.ko = pf(F1, length(ko), (length(wt)-3), lower.tail = FALSE, log.p = FALSE)
    }
    
  }
  
  #print(pval.wt.ko)
  return(c(pval.ko=pval.wt.ko, nb.ko=nb.ko, corr.ko =corr.ko))
}

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


model.sel.wt.ko.allModel = function(data.wt.ko, time.wt=c(0:15)*3, time.ko=c(0,6,12,18), period=24)
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
    prob.wt.ko = c(NA, NA, NA, NA, NA);
  }else{
    ### fitting wt.ko pool with same rhythmic parameters
    wt.ko = c(wt, ko)
    time.wt.ko = c(time.wt, time.ko)
    c=cos(2*pi*time.wt.ko/period)
    s=sin(2*pi*time.wt.ko/period)
    fit = lm(wt.ko~c+s)
    rss = sum(fit$residuals^2)
    
    ###fitting wt and ko with same flat parameter
    rss5 = sum((wt.ko-mean(wt.ko))^2)
    
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
    
    ###fitting wt and ko with same flat parameter
    rss4 = sum((wt-mean(wt))^2)
    
    ## model 1: rhythmic with same parameters
    n = length(wt.ko)
    BIC1 = n*log(rss/n) + 3*log(n);
    BIC2 = n*log((rss1+rss2)/n) + 6*log(n);
    BIC3 = n*log((rss1+rss3)/n) + 4*log(n);
    BIC4 = n*log((rss4+rss2)/n) + 4*log(n);
    BIC5 = n*log(rss5/n) + 1*log(n);
    
    BIC = c(BIC1, BIC2, BIC3, BIC4, BIC5)
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    
    prob.wt.ko = prob.model
  }
  
  #print(pval.wt.ko)
  names(prob.wt.ko) = c('prob.M1', 'prob.M2', 'prob.M3', 'prob.M4', 'prob.M5')
  return(c(prob.wt.ko, nb.ko=nb.ko, corr.ko =corr.ko))
}

model.sel.allModel = function(data1, t1=c(0:15)*3, data2, t2=c(0:15)*3, period=24)
{
  # period = 24; data1 = nuclear[3, c(1:16)]; data2 = nuclear[4, c(1:16)]; t1 = c(0:15)*3; t2 = t1;
  data1 = as.numeric(data1)
  data2 = as.numeric(data2)
  d1 = (data1 - mean(data1, na.rm = TRUE))/sd(data1, na.rm = TRUE);
  d2 = (data2 - mean(data2, na.rm = TRUE))/sd(data2, na.rm = TRUE);
  d1 = d1[which(!is.na(d1)==TRUE)];t1 = t1[which(!is.na(d1)==TRUE)];
  d2 = d2[which(!is.na(d2)==TRUE)];t2 = t2[which(!is.na(d2)==TRUE)];
  
  ### M1: fitting the pool of d1 and d2 with same rhythmic parameters
  d = c(d1, d2)
  t = c(t1, t2)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit = lm(d~c+s)
  rss = sum(fit$residuals^2)
  
  ### M2: fitting d1 and d2 with different rhythmic parameters
  c1=cos(2*pi*t1/period)
  s1=sin(2*pi*t1/period)
  fit1 = lm(d1~c1+s1)
  rss1 = sum(fit1$residuals^2)
  
  c2=cos(2*pi*t2/period)
  s2=sin(2*pi*t2/period)
  fit2 = lm(d2~c2+s2)
  rss2 = sum(fit2$residuals^2)
  
  ## model 1: rhythmic with same parameters
  n = length(d);
  BIC1 = n*log(rss/n) + 2*log(n);
  BIC2 = n*log((rss1+rss2)/n) + 4*log(n);
  #BIC1 = n*log(rss/n) + 2*2;
  #BIC2 = n*log((rss1+rss2)/n) + 4*2;
  
  BIC = c(BIC1, BIC2)
  bic = BIC-min(BIC)
  prob.model = exp(-0.5*bic)
  prob.model = prob.model/sum(prob.model)
  #prob.model = c(prob.model, which(prob.model==max(prob.model)))
  #print(pval.wt.ko)
  names(prob.model) = c('prob.BIC.M1', 'prob.BIC.M2')
  
  #plot(t1, d1, type='b', col='blue', ylim=range(c(d1, d2)))
  #points(t2, d2, type='b', col='black')
  
  return(prob.model)
  
}


########################################################
########################################################
# Section : protein complex analysis
# 
########################################################
########################################################
mean.correlation.rhythmic = function(kk, index)
{
  correlation = c()
  for(k in kk)
  {	
    correl = c()
    v1 = as.numeric(nuclear[k, c(1:16)]); 
    index.cor = index[which(index!=k)]
    for(ii in index.cor) 
    {	
      #print(c(k, ii))
      v2 = as.numeric(nuclear[ii,c(1:16)])
      correl = c(correl, mutual.correlation(v1, v2))
    }
    correlation = c(correlation, mean(correl))
  }
  return(correlation)
  
}
mean.correlation.subunits = function(index)
{
  correlation = c()
  ii = c(1:(length(index)-1))
  for(i in ii)
  {
    jj = c(1:length(index))
    jj = jj[which(jj>i)]
    for(j in jj)
    {
      v1 = as.numeric(nuclear[index[i], c(1:16)]); 
      v2 = as.numeric(nuclear[index[j], c(1:16)])
      correlation = c(correlation, mutual.correlation(v1, v2))
    }
  }
  return(mean(correlation))
}

##########################################
# curate protein complex annotation
##########################################
Processing.CORUM.all.Complexes.add.manually = function()
{
  annotation = read.csv('Annotations/allComplexes_CORUM.csv', sep=';')
  #ii = which(annotation$organism=='Mouse'|annotation$organism=='Human'|annotation$organism=='Rat')
  annot =annotation[,c(1,2,3,4,5)]
  
  mapping = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Protein_Complexes/refGene2uniprotID.txt', header=FALSE, sep='\t')
  uniprot = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Uniprot_proteins_ID.csv',header=TRUE)
  uniprot.mouse = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Protein_Complexes/uniprot-taxonomy-mouse.tab', header=TRUE)
  uniprot.human = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Protein_Complexes/uniprot-taxonomy-human.tab', header=TRUE)
  uniprot.rat = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Protein_Complexes/uniprot-taxonomy-rat.tab', header=TRUE)
  uniprot = rbind(uniprot.mouse, uniprot.human, uniprot.rat)
  
  annot.corum = annot;
  
  subunits = c() 
  nb.subunits = c()
  for(n in 1:nrow(annot.corum))
  {
    cat(n, '\n');
    test = annot.corum[n,5]
    test = gsub("[()]","",test) 
    test = unlist(strsplit(as.character(test),","))
    ii = match(test, uniprot[,1])
    if(length(which(!is.na(ii)==TRUE))!=length(test)) print(n)
    ii = ii[which(!is.na(ii)==TRUE)]
    nb.subunits = c(nb.subunits, length(ii))
    gene = c()
    for(kk in ii) 
    {
      gg = sapply(unlist(strsplit(as.character(uniprot$Gene.names[kk]), ' ')), first.upper, USE.NAMES = FALSE);
      mm = match(gg, microsomal.names[,3])
      if(length(which(!is.na(mm)==TRUE))>0) gene = c(gene, gg[which(!is.na(mm)==TRUE)][1])
      else gene = c(gene, gg[1])
    }
    subunits = c(subunits, paste(gene, sep='', collapse=','))
  }
  
  annot.corum = data.frame(annot.corum, nb.subunits, subunits, stringsAsFactors = FALSE)
  
  #### Mannually Add complexes
  ### Bmal1-Clock complex
  pc.m = c('C1', 'Bmal1-Clock Heterodimer', rep(NA, 5), 'Arntl,Clock')
  ### Per complex
  pc.m = rbind(pc.m, c('C2', 'PER complex (main)', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2'))
  pc.m = rbind(pc.m, c('C3', 'Bmal-Clock-Pers complex 1', rep(NA, 5), 'Arntl,Clock,Per1,Per2,Per3,Cry1,Cry2'))
  pc.m = rbind(pc.m, c('C4', 'Clock-Bmal1 complex 1', rep(NA, 5), 'Clock,Arntl,Gnb2l1,Prkaca,Thrap3')) ## Robles et al. 2010 and Lande-Diner et al. 2013
  pc.m = rbind(pc.m, c('C5', 'PER complex 1', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5')) ##Brown et al.2005
  pc.m = rbind(pc.m, c('C6', 'PER complex 2', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5,Csnk1d,Csnk1e,Sfpq,Sin3a,Sin3b,Hdac1,Hdac2,Sap18,Sap30,Rbbp4,Rbbp7')) ##Duong et al. 2011
  pc.m = rbind(pc.m, c('C7', 'PER complex 3', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5,Csnk1d,Csnk1e,Ddx5,Dhx9,Setx')) ## Padmanabhan et al. 2012
  pc.m = rbind(pc.m, c('C8', 'PER complex 4', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5,Csnk1d,Csnk1e,Suv39h1,Suv39h2,Cbx3,Trim28')) ## Duong et al. 2014
  
  pc.m = rbind(pc.m, c('C20', 'Bmal1-Clock complex 2', rep(NA, 5), 'Arntl,Clock,Mta2,Chd4')) ## Kim et al. 2014
  pc.m = rbind(pc.m, c('C21', 'PER complex 5', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Hdac1,Hdac2,Mbd2,Gatad2a,Rbbp4')) ## Kim et al. 2014
  pc.m = rbind(pc.m, c('C22', 'PER complex 6', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5,Csnk1d,Csnk1e,Suv39h1,Suv39h2,Cbx3,Trim28')) ## Kim et al. 2014
  
  ### Polymerase complexes
  pc.m = rbind(pc.m, c('C9', 'RNA Polymerase I complex (parts)', rep(NA, 5), 'Polr1a,Polr1b,Polr1c,Polr1d,Polr1e,Polr2e')) ## from Fred
  pc.m = rbind(pc.m, c('C10', 'RNA Polymerase I complex (whole)', rep(NA, 5), 'Polr1a,Polr1b,Polr1c,Polr1d,Polr1e,Polr2e,Polr2f,Polr2h,Polr2l,Polr2k,Polr1e'))
  pc.m = rbind(pc.m, c('C11', 'RNA Polymerase II complex (whole)', rep(NA, 5), 'Polr2a,Polr2b,Polr2c,Polr2j,Polr2l,Polr2e,Polr2f,Polr2h,Polr2l,Polr2k,Polr2d,Polr2g,Gtf2f1,Gtf2f2'))
  pc.m = rbind(pc.m, c('C12', 'RNA Polymerase III complex (whole)', rep(NA, 5), 'Polr3a,Polr3b,Polr1c,Polr1d,Polr3d,Polr2e,Polr2f,Polr2h,Polr2l,Polr2k,Polr3h,Polr3k,Polr3c,Polr3f,Polr3g,Polr3l,Polr3e'))
  
  ### add methyltransferase complexes
  pc.m = rbind(pc.m, c('C13', 'H3k4 methyltransferase complex (Set1a)', rep(NA, 5), 'Setd1a,Ash2l,Rbbp5,Wdr5,Dpy30,Cxxc1,Wdr82,Hcfc1,Hcfc2'))
  pc.m = rbind(pc.m, c('C14', 'H3k4 methyltransferase complex (Set1b)', rep(NA, 5), 'Setd1b,Ash2l,Rbbp5,Wdr5,Dpy30,Cxxc1,Wdr82,Bod1,Bod1l,Hcfc1,Hcfc2'))
  pc.m = rbind(pc.m, c('C15', 'H3k4 methyltransferase complex (Mll1)', rep(NA, 5), 'Mll1,Ash2l,Rbbp5,Wdr5,Dpy30,Hcfc1,Hcfc2,Men1'))
  pc.m = rbind(pc.m, c('C16', 'H3k4 methyltransferase complex (Mll2)', rep(NA, 5), 'Mll2,Ash2l,Rbbp5,Wdr5,Dpy30,Hcfc1,Hcfc2,Men1,Psip1'))
  pc.m = rbind(pc.m, c('C17', 'H3k4 methyltransferase complex (Mll3)', rep(NA, 5), 'Mll3,Ash2l,Rbbp5,Wdr5,Dpy30,Ncoa6,Kdm6a,Paxpip1,Pagr1'))
  pc.m = rbind(pc.m, c('C18', 'H3k4 methyltransferase complex (Mll4)', rep(NA, 5), 'Mll4,Ash2l,Rbbp5,Wdr5,Dpy30,Ncoa6,Kdm6a,Paxpip1,Pagr1'))
  pc.m = rbind(pc.m, c('C19', 'H3k4 methyltransferase complex (Mll5)', rep(NA, 5), 'Mll5,Hcfc1,Ogt,Stk38,Ppp1ca,Ppp1cb,Ppp1cc,Actb'))
  
  pc.m = rbind(pc.m, c('C23', 'HIRA histone chaperone complex', rep(NA, 5), 'Hira,Ubn1,Cabin1,Asf1a'))
  pc.m = rbind(pc.m, c('C24', 'Rev-ErbA-alpha-Ncor1-Hdac3', rep(NA, 5), 'Nr1d1,Ncor1,Hdac3')) ## Yin and Lazar, 2005
  
  #pc.m = rbind(pc.m, c('C13', 'PER complex', rep(NA, 5), 'Per1,Per2,Per3,Cry1,Cry2,Nono,Wdr5,Csnk1d,Csnk1e,Sfpq'))
  pc.m = pc.m[, c(1:6, 8)]
  colnames = colnames(annot.corum)
  colnames(pc.m) = colnames
  pc.m = data.frame(pc.m, stringsAsFactors=FALSE)
  
  nb = c()
  for(n in 1:nrow(pc.m))
  {
    ss = pc.m$subunits[n]
    ss = unlist(strsplit(as.character(ss), ',')) 
    nb = c(nb, length(ss))
  }
  pc.m$nb.subunits = nb
  xx = rbind(pc.m, annot.corum)
  #colnames(xx) = colnames
  
  annot.corum = xx
  save(annot.corum, file='Rdata/Annotation_Corum_all_with_Manual.Rdata')
  #write.table(uniprot, file='uniprot-GeneIDs-mapping-human-mouse-rat_v2.txt',quote=FALSE, sep='\t',col.names=TRUE, row.names=FALSE)
  #write.table(annot.corum, file='Annotation_complexes_corum_v2.txt',quote=FALSE, sep='\t',col.names=TRUE, row.names=FALSE) 
  return(annot.corum)
}

statistics.complexes = function(vect)
{
	cutoff.rhy = 0.01
	cutoff.sim = 0.80
	cutoff.nb.timepoints = 8
	
	stat = c()
	
	index = vect[which(names(vect)=='index.detected')];
	index = as.numeric(unlist(strsplit(as.character(index), ',')))
#index = index[which(nuclear$nb.timepoints[index]==16)];
	index = index[which(nuclear$nb.timepoints[index]>=cutoff.nb.timepoints)]
	stat = c(stat, length(index))
	
	if(length(index)<2){
		stat = c(stat, rep(NA, 9))
	}else{			
		
		stat = c(stat, paste(nuclear$pval[index], sep='', collapse=','), paste(nuclear$qv[index], sep='', collapse=','), paste(nuclear$phase[index], sep='', collapse=','))
		
		kk = index[which(nuclear$pval[index]<cutoff.rhy)]
		stat = c(stat, length(kk), length(kk)/length(index))
		
		if(length(kk)>0){
			
			stat = c(stat, paste(nuclear$phase[kk], sep='', collapse=','));
### phase differences for rhythmic subunits
			if(length(kk)==1){
				stat = c(stat, NA);
			}else{
				diff.phase = phase.difference(kk)
				stat = c(stat, paste(diff.phase, sep='', collapse=','))
			}
### correlation of rhythmic subunits with other subunits
			stat = c(stat, paste(mean.correlation.rhythmic(kk, index), sep='', collapse=','));
			
		}else{
			stat = c(stat, NA, NA, NA);
		}
### mean.correlation among subunits
		stat = c(stat, mean.correlation.subunits(index));
		
	}
	
	names(stat) = c('nb.subunits.quantified','pval.subunits', 'qv.subunits', 'phase.subunits','nb.rhythmic.subunits', 'percentage.rhythmic.subunits','phase.rhythmic.subunits', 'phase.diff.rhythmic.subunits',  'correlation.rhythmic.subunits', 'mean.correlation.subunits');
	return(stat);
}

statistics.complexes.all.fitting = function(vect)
{
	cutoff.rhy = 0.01
	cutoff.sim = 0.80
	cutoff.nb.timepoints = 8
	
	stadardization.nona = function(x)
	{	
		xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
		return(xx)
	}
	
	stat = c()
	
	index = vect[which(names(vect)=='index.detected')];
	index = as.numeric(unlist(strsplit(as.character(index), ',')))
  #index = index[which(nuclear$nb.timepoints[index]==16)];
	index = index[which(nuclear$nb.timepoints[index]>=cutoff.nb.timepoints)]
	stat = c(stat, length(index))
	
	if(length(index)<2){
		stat = c(stat, rep(NA, 6))
	}else{			
		test = c()
		for(iidex in index) test = c(test, stadardization.nona(as.numeric(nuclear[iidex, c(1:16)])))
		stat = c(stat, f24_R2_alt2(test, t=c(0:(length(test)-1)*3)));
		
	}
	
	names(stat)[1] = 'nb.subunits.quantified'
	return(stat);
}

##########################################
# protein complex plots
##########################################
Plots.complexes.all.fitting = function(kk, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
	pdf(pdfname, width = 14, height = 12)
	stadardization.nona = function(x)
	{	
#xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
		xx = x-mean(x[which(!is.na(x)==TRUE)]);
		return(xx)
	}
	
	ylim = c(-2.5, 2.5)
	xlim = c(0,48)
	
	for(n in kk)
	{
		
		index = res$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(res$subunits.detected[n]), ','))
		
		index = index[which(nuclear$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, stadardization.nona)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			ttest = as.vector(test)
			ttest =ttest[which(!is.na(test))]
			ylim = range(ttest)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(res[n,2],'\n', 'pval.all=', signif(as.numeric(res$pval[n]), d=3),																																		
																													 ', qval.all=', signif(as.numeric(res$qval[n]),d=3), 
																													 ', phase.all=', signif(as.numeric(res$phase[n]), d=3)))
			abline(h=0, col='gray',lwd=4.0)
			
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=1.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}

Plots.complexes.variance_cluster = function(kk, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
	pdf(pdfname, width = 14, height = 12)
	stadardization.nona = function(x)
	{	
		#xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
		xx = x-mean(x[which(!is.na(x)==TRUE)]);
		return(xx)
	}
	
	ylim = c(-2.5, 2.5)
	xlim = c(0,48)
	
	for(n in kk)
	{
		
		index = res$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(res$subunits.detected[n]), ','))
		
		index = index[which(nuclear$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, stadardization.nona)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			ttest = as.vector(test)
			ttest =ttest[which(!is.na(test))]
			ylim = range(ttest)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(res[n,2],'\n', 'pval.all=', signif(as.numeric(res$pval[n]), d=3),																																		
																													 ', qval.all=', signif(as.numeric(res$qval[n]),d=3), 
																													 ', phase.all=', signif(as.numeric(res$phase[n]), d=3)))
			abline(h=0, col='gray',lwd=4.0)
			
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=1.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}



Plots.complexes = function(kk, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
	pdf(pdfname, width = 14, height = 12)
		
	ylim = c(-1.5, 1.5)
	xlim = c(0,48)
	
	for(n in kk)
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		index = index[which(nuclear$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, standadization.nona.impute)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			ylim = range(test, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],',\n', 'Coverage=', signif(as.numeric(annot$percent.detected[n]),d=2)*100, '%'))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}

Plots.complexes.v2 = function(annot, kk, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
	pdf(pdfname, width = 14, height = 12)
	
	ylim = c(-1.5, 1.5)
	xlim = c(0,48)
	
	for(n in kk)
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		index = index[which(nuclear$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, standadization.nona.impute)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			ylim = range(test, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],',\n', 'Coverage=', signif(as.numeric(annot$percent.detected[n]),d=2)*100, '%'))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}

Plots.complexes.v3 = function(annot, nn, pdfname)
{
	table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v2.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
	
	kk = which(!is.na(table.sx$TFs)==TRUE)
	tfs = table.sx[kk,2]
	kk = which(!is.na(table.sx$Transcription.Cofactors)==TRUE)
	cofactors = table.sx[kk,2]
	kk = which(!is.na(table.sx$RNA.Processing)==TRUE)
	rna.processing.protein = table.sx[kk,2]
	
	pdf(pdfname, width = 14, height = 12)
	
	ylim = c(-1.5, 1.5)
	xlim = c(0,48)
	
	for(n in nn)
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		ii = which(nuclear$nb.timepoints[index]>=8)
		index = index[ii]
		genes = genes[ii]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, standadize.nona)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=2)
			amps = signif(nuclear$amp[index],d=2)
			
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			amps = amps[o1]
			index = index[o1]
			genes = genes[o1]
			ylim = range(test, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],',\n', 'Coverage=', signif(as.numeric(annot$percent.detected[n]),d=2)*100, '%', 
																													 ', pval.svd = ', signif(as.numeric(annot$pval.svd[n]), d=2), ', qv.rhythmic = ', signif(as.numeric(annot$qv.p1[n]), d=2), sep=''))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], cofactors)
				if(!is.na(mm)) {
					text.col = 'lightslateblue'
				}else{
					if(!is.na(match(genes[i], tfs))){
					text.col='magenta'
					}else{
						if(!is.na(match(genes[i], rna.processing.protein))){
							text.col = 'orange'
							}else{
								text.col = 'black'
							}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', amp= ', amps[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}

Plots.complexes.subunits = function(annot, nn, prot.data, pdfname)
{
  #table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v2.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
  #load(file='Rdata/Table_microsomal_log2_L_H_names_mRNA_total_nuclear.Rdata')
  nuclear = prot.data
  
  pdf(pdfname, width = 14, height = 12)
  
  ylim = c(-1.5, 1.5)
  xlim = c(0,48)
  
  for(n in nn)
  {
    cat(n, '\n')
    index = annot$index.detected[n]
    index = as.numeric(unlist(strsplit(as.character(index), ',')))
    genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
    #genes = genes[which(nuclear$nb.timepoints[index]==16)]
    ii = which(nuclear$nb.timepoints[index]>=8)
    index = index[ii]
    genes = genes[ii]
    
    if(length(index)>1)
    {
      rainbow = rainbow(length(index),s = 0.85, v = 0.85)
      
      test = as.matrix(nuclear[index,c(1:16)])
      test =apply(test, 1, standadize.nona)
      test = t(test)
      pvals = signif(nuclear$pval.WT.24hRhythmicity[index],d=2)
      phases = signif(nuclear$phase.WT.24hRhythmicity[index],d=2)
      amps = signif(nuclear$amp.WT.24hRhythmicity[index],d=2)
      
      o1 = order(pvals)
      pvals = pvals[o1]
      test = test[o1,]
      phases = phases[o1]
      amps = amps[o1]
      index = index[o1]
      genes = genes[o1]
      ylim = range(test, na.rm=TRUE)
      
      plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,
           main=paste(annot[n,2],',\n', 
                      #'Coverage=',  annot$percent.detected[n], 
                      'pval.svd = ', 
                      signif(as.numeric(annot$pval.svd[n]), d=2), 
                      ', pval.rhythmic = ', 
                      signif(as.numeric(annot$pval.p1[n]), d=2), sep='')
           )
      abline(h=0, col='gray',lwd=4.0)
      
      for(i in 1:nrow(test))
      {
        i.rel = i; 
        text.col = 'black';
        
        type.plot='b';
        points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, 
               col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
        #points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
        legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), 
               legend = paste(genes[i],', pval= ', pvals[i], ', amp= ', amps[i],sep=''),
               text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
      }
      
    }
  }
  dev.off()
  
}

Plots.complexes.details = function(annot, nn, pdfname)
{
	#regulator.all = c('tfs.curated', 'cofactors', 'chromatin.remodellers',  
  #'rna.processing.proteins','splicesome', 'kinases', 'phosphatases')
	#load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_kinases_phosphatases_used.Rdata')
	#tfs = tfs.curated
	#cofactors = unique(c(cofactors, chromatin.remodellers))
	
	table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v2.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
	
	kk = which(!is.na(table.sx$TFs)==TRUE)
	tfs = table.sx[kk,2]
	kk = which(!is.na(table.sx$Transcription.Cofactors)==TRUE)
	cofactors = table.sx[kk,2]
	kk = which(!is.na(table.sx$RNA.Processing)==TRUE)
	rna.processing.protein = table.sx[kk,2]

	
	#tfs = tfs[,1]
	pdf(pdfname, width = 16, height = 12)
	par(mfrow=c(2,2))
	
	ylim = c(-1.5, 1.5)
	xlim = c(0,50)
	
	for(n in nn)
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		ii = which(nuclear$nb.timepoints[index]>=8)
		index = index[ii]
		genes = genes[ii]
		bg.cutoff = annot$bg.cutoff[n]
		
		if(length(index)>1)
		{
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, standadize.nona)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=2)
			amps = signif(nuclear$amp[index],d=2)
			
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			amps = amps[o1]
			index = index[o1]
			genes = genes[o1]
			
			#### Panel A: profiles of subunits
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			ylim = range(test, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'a.u', xlim = xlim, ylim = ylim,main=paste(annot[n,2],',\n', sep=''))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], cofactors)
				if(!is.na(mm)) {
					text.col = 'lightslateblue'
				}else{if(!is.na(match(genes[i], tfs))){
						text.col='magenta'
					}else{
						if(!is.na(match(genes[i], rna.processing.protein))){
						text.col = 'orange'
						}else{
							text.col = 'black'
						}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.85*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/10,(1.0/nrow(test)))), legend = paste(genes[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
			
			#### Panel B: amplitudes and p values of subunits
			ppvals = -log10(pvals)
			plot(ppvals, amps, type='p', pch=19, col=rainbow, cex=1.5, xlim=c(0, max(ppvals)), ylim=c(0, max(amps)), xlab='-log10(pvals)')
			abline(v=2, col='darkgray', lwd=2.0);
			abline(h=0.1, col='darkgray', lwd=2.0);
			##### Panel C: components of SVD
			test = matrix(data = NA, nrow = length(index), ncol = 16)
			for(kk in 1:length(index))
			test[kk, ] = standadization.nona.impute(as.numeric(nuclear[index[kk], c(1:16)]))
			#rownames(test) = subunits
			###SVD
			ss = svd(test)
			L = length(ss$d)
			pl = ss$d^2/sum(ss$d^2)
						
			#### first 2 components, first change the first two conlums to positive
			u1 = ss$u[,1]
			u2 = ss$u[,2]
			v1 = ss$v[,1]
			v2 = ss$v[,2]
			nb.positive = length(which(u1>=0))
			nb.negative = length(which(u1<0))
			if(nb.positive<nb.negative){u1 = -u1; v1 = -v1;}
			
			nb.positive = length(which(u2>=0))
			nb.negative = length(which(u2<0))
			if(nb.positive<nb.negative){u2 = -u2; v2 = -v2;}
			
			res1 = pl[1]*length(which(u1>=0))/length(u1);
			stat1 = f24_R2_alt2(v1, t=c(0:15)*3)
			
			plot(c(1:length(pl)), pl, type='b', pch=16, lwd=2.0,cex=1.0, ylab='% variance explained by each component', xlab='Index of components', ylim=c(0,1),main=NA)
			abline(h=bg.cutoff, col='darkgray',lwd=2.5)
			
			plot(c(0:15)*3, apply(ss$u, 2, mean)[1]*ss$v[,1], type='b',lwd=2.0,col='darkblue',ylab='1st component',main=paste('amp=', signif(stat1[3],d=2), ', phase= ', signif(stat1[5], d=2), ', pval=', signif(stat1[6], d=3)))
			abline(h=0, col='darkgray',lwd=2.0)
						
		}
	}
	dev.off()
	
}






Plots.complexes.wt.ko = function(aa, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
		
	ylim = c(-1.5, 1.5)
	xlim = c(0,48)
	
	for(n in 1:nrow(aa))
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		index = index[which(nuclear.all$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			
			
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear.all[index,c(1:16)])
			test.m = as.matrix(nuclear.all[index, c(34:45)])
			test.m1 = test.m[,c(1:4)]
			test.m2 = test.m[,c(5:8)]
			test.m3 = test.m[,c(9:12)]
			#test =apply(test, 1, centering.nona)
			#test = t(test)
			#test.m1 = t(apply(test.m1, 1, centering.nona))
			#test.m2 = t(apply(test.m2, 1, centering.nona))
			#test.m3 = t(apply(test.m3, 1, centering.nona))
			pvals = signif(nuclear.all$pval[index],d=2)
			phases = signif(nuclear.all$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			
			pdfname = paste('myplots/PCs_circadian_clock_WT_KO_', aa[n,2], '.pdf')
			pdf(pdfname, width = 12, height = 10)
			par(mfcol = c(2,2))
				
			## WI
			ylim = range(cbind(test, test.m1, test.m2, test.m3), na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'centered abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],',\n', 'Coverage=', signif(as.numeric(annot$percent.detected[n]),d=2)*100, '%'))
			abline(h=0, col='gray',lwd=2.0)
			cex = 1.5;
			lwd= 2.0;
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=cex, lwd=lwd, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.6*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/10,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			
			}
			## Cry ko
			#ylim = range(test.m1, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'centered abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],', CRY KO'))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				type.plot='b';
				points(c(0:3)*6, (test.m1[i,]), type=type.plot, cex=cex, lwd=lwd, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				legend(x = 0.6*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/10,(1.0/nrow(test)))), legend = paste(genes[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
			## Bmal WT
			#ylim = range(test.m2, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'centered abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],', Bmal WT'))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				type.plot='b';
				points(c(0:3)*6, (test.m2[i,]), type=type.plot, cex=cex, lwd=lwd, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				legend(x = 0.6*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/10,(1.0/nrow(test)))), legend = paste(genes[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
			## Bmal KO
			#ylim = range(test.m3, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'centered abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2],', Bmal KO'))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{	
				i.rel = i; 
				type.plot='b';
				points(c(0:3)*6, (test.m3[i,]), type=type.plot, cex=cex, lwd=lwd, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				legend(x = 0.6*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/10,(1.0/nrow(test)))), legend = paste(genes[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
			
			dev.off()	
		}
	}
	
	
}


stadardization.nona = function(x)
{	
  xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
  return(xx)
}

statistics.complexes.all.fitting.2 = function(vect)
{
	cutoff.rhy = 0.01
	cutoff.sim = 0.80
	cutoff.nb.timepoints = 8
	
	
	
	stat = c()
	
	index = vect[which(names(vect)=='index.detected')];
	index = as.numeric(unlist(strsplit(as.character(index), ',')))
	#index = index[which(nuclear$nb.timepoints[index]==16)];
	index = index[which(nuclear$nb.timepoints[index]>=cutoff.nb.timepoints)]
	stat = c(stat, length(index))
	
	if(length(index)<2){## the cases where less than 2 subunits are detected
		stat = c(stat, rep(NA, 8))
	}else{ ## the cases where more than 2 subunits are detected
		test = matrix(data = NA, nrow = length(index), ncol = 16)
		for(kk in 1:length(index)) test[kk, ] = stadardization.nona(as.numeric(nuclear[index[kk], c(1:16)]))
		
		var.pc = c()
		mean.pc = c()
		for(kk in 1:16)
		{
			ttest = test[,kk]
			ttest = ttest[which(!is.na(ttest))]
			if(length(ttest)<1){var.pc = c(var.pc, NA);mean.pc = c(mean.pc, NA);	}
			if(length(ttest)==1){var.pc = c(var.pc, NA);mean.pc = c(mean.pc, mean(ttest));}
			if(length(ttest)>1){var.pc = c(var.pc, var(ttest));mean.pc = c(mean.pc, mean(ttest));}
		}
		var.pc = var.pc[which(!is.na(var.pc))]
		if(length(var.pc)>=1) {
			stat = c(stat, f24_R2_alt2(mean.pc, t=c(0:15)*3), mean(var.pc), max(var.pc));
		}else{
			stat = c(stat, f24_R2_alt2(mean.pc, t=c(0:15)*3), NA, NA);
		}
		
	}
	return(stat);
}

#####
##### Function to quantify the subunits forming complexe and when they form complex 
#####
centering.nona = function(x)
{	
	xx = x-mean(x[which(!is.na(x)==TRUE)])
	return(xx)
}
mean.nona = function(x)
{	
  kk = which(!is.na(x))
  if(length(kk)>0) 
  {
    xx = mean(x[which(!is.na(x)==TRUE)]);
  }else{
    xx = NA;
  }
	return(xx)
}
sme.nona = function(x)
{
    kk = which(!is.na(x))
    if(length(kk)>1)
    {
        return(sd(x[kk])/sqrt(length(kk)))
    }else{
        return(NA)
    }
}

standadize.nona = function(xx)
{	
	kk = which(!is.na(xx)==TRUE);
	xx = (xx-mean(xx[kk]))/sd(xx[kk])
	return(xx)
}


compute.sd = function(s1, s2)
{
	kk = intersect(which(!is.na(s1)==TRUE), which(!is.na(s2)==TRUE))
	delta.s = s1[kk] - s2[kk]
	return(sd(delta.s))
}
paste.vector = function(vect)
{
	return(paste(vect, collapse=','))
}

centering.nona.impute = function(x)
{	
	kk = which(is.na(x)==TRUE);
	xx = x-mean(x[which(!is.na(x)==TRUE)])
	if(length(kk)>0)
	{
		stat = f24_R2_alt2(xx, t=c(0:15)*3)
		xx[kk] = stat[2] + stat[3]/2*cos(2*pi/24*((kk-1)*3-stat[5]))  
	}
	return(xx)
}
standadization.nona.impute = function(x)
{	
	kk = which(is.na(x)==TRUE);
	xx = x
	if(length(kk)>0)
	{
		stat = f24_R2_alt2(xx, t=c(0:15)*3)
		xx[kk] = stat[2] + stat[3]/2*cos(2*pi/24*((kk-1)*3-stat[5]))  
	}
	xx = (xx-mean(xx))/sd(xx)
	return(xx)
}



statistics.complexes.svd = function(annot, prot.data, res, pdfname='SVD_PC_plot_example.pdf', 
                                    index.data = c(1:16), 
                                    cutoff.nb.timepoints = 8, TEST = FALSE)
{	
	colnames(res) = c('nb.subunits.quantified', 'd.entropy', 'svd.1st.component.p1','svd.2nd.compoent.p2', 
	                  'amp.p1', 'phase.p1','pval.p1', 'amp.p2', 'phase.p2','pval.p2')
	## ii = grep('RNA polymerase II holoenzyme complex', annot[,2])[1];vect = annot[ii,]
	if(TEST){annot = annot[c(1:20), ]; }
	
	pdf(pdfname, width=14, height=12)
	
	for(n in 1:nrow(annot))
	{
    # n = 1
		cat(n, '\n');
		vect = annot[n,]
		index = vect[which(names(vect)=='index.detected')];
		subunits = vect[which(names(vect)=='subunits.detected')];
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		subunits = unlist(strsplit(as.character(subunits), ','))
		
		ii.sels = which(prot.data$nb.timepoints[index]>=cutoff.nb.timepoints)
		subunits = subunits[ii.sels]
		index = index[ii.sels]
		
		res[n,1] = length(index)
		
		if(length(index)>=2) {
			test = matrix(data = NA, nrow = length(index), ncol = 16)
			for(kk in 1:length(index))
			test[kk, ] = standadization.nona.impute(as.numeric(prot.data[index[kk], index.data]))
			rownames(test) = subunits
			
			# SVD
			ss = svd(test)
			L = length(ss$d)
			pl = ss$d^2/sum(ss$d^2)
			dd = -1/log(L)*sum(pl*log(pl)) 
			res[n, 2] = dd;
						
			#### first 2 components, first change the first two conlums to positive
			u1 = ss$u[,1]
			u2 = ss$u[,2]
			v1 = ss$v[,1]
			v2 = ss$v[,2]
			nb.positive = length(which(u1>=0))
			nb.negative = length(which(u1<0))
			if(nb.positive<nb.negative){u1 = -u1; v1 = -v1;}
			
			nb.positive = length(which(u2>=0))
			nb.negative = length(which(u2<0))
			if(nb.positive<nb.negative){u2 = -u2; v2 = -v2;}
			
			res[n, 3] = pl[1]*length(which(u1>=0))/length(u1);
			res[n, 4] = pl[2]*length(which(u2>=0))/length(u2);
			
			stat1 = f24_R2_alt2(v1, t=c(0:15)*3)
			stat2 = f24_R2_alt2(v2, t=c(0:15)*3)
			res[n, 5] = stat1[3]
			res[n, 6] = stat1[5]
			res[n, 7] = stat1[6]
			res[n, 8] = stat2[3]
			res[n, 9] = stat2[5]
			res[n, 10] = stat2[6]
			#print(pl)
			#print(dd)
			#par(mfcol = c(1,2))
			#matplot(t(X), type='b')
			#matplot(s$v[,c(1:3)], type='b')
			par(mfcol = c(2,2))
			matplot(c(0:15)*3, t(test), type='b',lwd=1.5, ylab='Abundance of subunits', main=annot[n,2])
			abline(h=0, col='gray',lwd=2.0)
			
			plot(c(1:length(pl)), pl, type='b', pch=16, lwd=2.0,cex=1.0, 
			     ylab='% variance explained by each component', 
			     xlab='Index of components', ylim=c(0,1), 
			     main=paste('d= ', signif(dd, d=2),sep=''))
			abline(v=(which(pl>0.7/length(pl))[length(which(pl>0.7/length(pl)))]+0.2), col='red',lwd=2.5)
			
			plot(c(0:15)*3, apply(ss$u, 2, mean)[1]*ss$v[,1], type='b',lwd=2.0,col='darkblue',
			     ylab='1st component',main=paste('amp=', signif(stat1[3],d=2), ', phase= ', signif(stat1[5], d=2), 
			                                     ', pval=', signif(stat1[6], d=3)))
			abline(h=0, col='gray',lwd=2.0)
			plot(c(0:15)*3, apply(ss$u, 2, mean)[2]*ss$v[,2], type='b',lwd=2.0, col='darkgreen', ylab='2nd component', 
			     main=paste('amp=', signif(stat2[3],d=2), ', phase= ', signif(stat2[5], d=2), ', pval=', signif(stat2[6], d=3)))
			abline(h=0, col='gray',lwd=2.0)
		}
	}
	
	dev.off()
	
	return(res)
}


##########################################
# clean the PC redundancy 
##########################################
reduce.redundancy.protein.complex = function(res)
{
  index.detected.clean = c()
  
  for(n in 1:nrow(res))
  {
    test = res$index.detected[n]
    test = unlist(strsplit(as.character(test), ','))
    test = as.numeric(test)
    test = unique(test)
    test = test[order(test)]
    index.detected.clean = c(index.detected.clean, paste(test, sep='', collapse=','))
  }
  
  res$index.detected.clean = index.detected.clean
  test = unique(res$index.detected.clean)
  index = c()
  id = c()
  names = c()
  species = c()
  percent.detected = c()
  
  for(n in 1:length(test))
  {
    jj = which(res$index.detected.clean==test[n])
    if(length(jj)==1)
    {
      index = c(index,jj)
      
      id = c(id, as.character(res[jj,1]))
      names = c(names, as.character(res[jj,2]))
      species = c(species, as.character(res[jj,4]))
      percent.detected = c(percent.detected, res$percent.detected[jj])
    }
    if(length(jj)>1) {
      print(n); 
      print(res[jj,c(1:2,4)]);
      print('......');
      
      index = c(index,jj[1])
      
      id = c(id, paste(res[jj,1], sep='', collapse=','))
      names = c(names, paste(res[jj,2], sep='', collapse=','))
      species = c(species, paste(res[jj,4],sep='', collapse=','))
      percent.detected = c(percent.detected, paste(res$percent.detected[jj], sep='', collapse=','))
    }
  }
  xx = data.frame(id, names, species, res[index,-c(1:4)])
  colnames(xx)[c(1:3)] = colnames(res)[c(1,2,4)]
  xx$percent.detected = percent.detected
  res = xx
  
  return(res)
  
}


Plots.complexes.svd = function(kk, pdfname)
{
	load(file='Rdata/Annotation_TFs_cofactors_chromatin_regulators_rna_processing_splicesome_used.Rdata')
	tfs = tfs[,1]
	pdf(pdfname, width = 14, height = 12)
	
	ylim = c(-1.5, 1.5)
	xlim = c(0,48)
	
	for(n in kk)
	{
		index = annot$index.detected[n]
		index = as.numeric(unlist(strsplit(as.character(index), ',')))
		genes = unlist(strsplit(as.character(annot$subunits.detected[n]), ','))
		#genes = genes[which(nuclear$nb.timepoints[index]==16)]
		index = index[which(nuclear$nb.timepoints[index]>=8)]
		
		if(length(index)>1)
		{
			rainbow = rainbow(length(index),s = 0.85, v = 0.85)
			test = as.matrix(nuclear[index,c(1:16)])
			test =apply(test, 1, standadization.nona.impute)
			test = t(test)
			pvals = signif(nuclear$pval[index],d=2)
			phases = signif(nuclear$phase[index],d=3)
			o1 = order(pvals)
			pvals = pvals[o1]
			test = test[o1,]
			phases = phases[o1]
			index = index[o1]
			genes = genes[o1]
			ylim = range(test, na.rm=TRUE)
			plot(1,1,type = 'n', xlab = 'ZT[h]', ylab = 'standadized abundance', xlim = xlim, ylim = ylim,main=paste(annot[n,2], ', Coverage=', signif(as.numeric(annot$percent.detected[n]),d=2)*100, '%, \ndd.entropy= ', signif(as.numeric(annot$d.entropy[n]), d=2), ', \np1= ', signif(as.numeric(annot$p1[n])*100, d=2), '%; p2= ', signif(as.numeric(annot$p2[n]), d=2)*100, '%, pval.p1= ',	signif(as.numeric(annot$pval.p1[n]), d=2), ', pval.p2= ', signif(as.numeric(annot$pval.p2[n]),d=2)))
			abline(h=0, col='gray',lwd=4.0)
			
			for(i in 1:nrow(test))
			{
				i.rel = i; 
				mm = match(genes[i], tfs)
				if(!is.na(mm)) {
					text.col = 'red'
				}else{
					jj = match(genes[i], cofactor)
					if(!is.na(jj)){
						text.col = 'blue'
					}else{
						if(!is.na(match(genes[i], chromatin.regulator))){
							text.col='green'
						}else{
							if(!is.na(match(genes[i], rna.processing.protein))){
								text.col='orange'
							}else{
								text.col = 'black'
							}
						}
					}
				}
				type.plot='b';
				points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
				#points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
				legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', phase= ', phases[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
			}
		}
	}
	dev.off()
	
}



simulation.phospho.nuclear = function(time, s0, epsilon.s0, phase.s0, s1, epsilon.s1, phase.s1, gamma0, epsilon.gamma0, phase.gamma0, gamma1, epsilon.gamma1, phase.gamma1,
									mu.kinase, epsilon.kinase, phase.kinase, mu.phosphatase, epsilon.phosphatase, phase.phosphatase, rate.kinase, rate.phosphatase, kappa.kinase, kappa.phosphatase)
{
	cat('Simulation')
	
	p0 = rep(0, length(time))
	p1 = rep(0, length(time))
	t = time
	
	Tstable =  24*(ceiling(log(2000)/max(gamma0, gamma1, s0, s1)/24) + 5) ## burning time is two period
	t.res = 2; 
	if(length(t)!=1){t.res = (t[2]-t[1])}
	t.sup = seq(0, Tstable+max(t),by= t.res)
	par = c(s0, epsilon.s0, phase.s0, s1, epsilon.s1, phase.s1, gamma0, epsilon.gamma0, phase.gamma0, gamma1, epsilon.gamma1, phase.gamma1,
			mu.kinase, epsilon.kinase, phase.kinase, mu.phosphatase, epsilon.phosphatase, phase.phosphatase, rate.kinase, kappa.kinase, rate.phosphatase, kappa.phosphatase)
	
	soln = lsoda(y = c(2, 4), #init.conditions
				 times= t.sup, ## times
				 dydt,
				 par=par, ## parameter values
				 rtol = 1e-6, atol = 1e-6) 
	
	#soln[match(t+48, soln[,1]),2]
	i.last = nrow(soln); 
	i.keep = seq(i.last - length(t)+1, i.last,by = 1)
	#plot(soln[,1],soln[,2], type = 'b',col='red')
	p0 = soln[i.keep,2]; 
	p1 = soln[i.keep,3]; 
	#mean = mean(m);m = m/mean
	#cat(soln[,2],'\n')
	#cat(parametrization,'\n')
	
	return(rbind(non.phospho=p0, phospho=p1, total=(p0+p1), ratio=(p1/(p0+p1))))
		
}

dydt = function(t, y, par)
{
	w = 2*pi/24
	p0 = y[1]
	p1 = y[2]
	
	#### sources of non-phospho and phospho
	s0 = par[1]*(1+par[2]*cos(w*(t-par[3])))
	s1 = par[4]*(1+par[5]*cos(w*(t-par[6])))
	#### degradation of non-phospho and phospho
	gamma0 = par[7]*(1+par[8]*cos(w*(t-par[9])))
	gamma1 = par[10]*(1+par[11]*cos(w*(t-par[12])))
	#### expression levels of kinase and phoshatase
	kk0 = par[13]*(1+par[14]*cos(w*(t-par[15])))
	kk1 = par[16]*(1+par[17]*cos(w*(t-par[18])))
	#### activities of kinase and phosphatase
	k0 = par[19]*kk0/(kk0 + par[20])
	k1 = par[21]*kk1/(kk1 + par[22])
		
	dp0dt = s0 - k0*p0 + k1*p1 - gamma0*p0
	dp1dt = s1 + k0*p0 - k1*p1 - gamma1*p1
	
	list(c(dp0dt, dp1dt),NULL)
}

difference.circular = function(p2, p1, period=24)
{
	if(p2<p1) p2 = p2 + period
	p = p2 - p1
	
	if(p<0) p=p+period
	if(p>period) p=p-period
	
	return(p)
	
}
mean.circular = function(p1, p2, period = 24)
{
	pp = c(p1, p2)
	pp = pp/period*2*pi
	a = mean(cos(pp))
	b = mean(sin(pp))
	p =period/(2*pi)*atan2(b, a)
	if(p<0) p=p+period
	if(p>period) p=p-period
	
	return(p)
	
}

mean.err = function(x, period=24, interval=3)
{
    x = as.numeric(x)
    kk = which(!is.na(x)==TRUE)
    index = period/interval;
    nn = length(x)/(period/interval)
    test = c()
    for(mm in 1:nn) test = cbind(test, x[c(1:index)+index*(mm-1)])
    test = data.frame(test)
    
    mean = apply(test, 1, mean.nona)
    err = apply(test, 1, sme.nona)
    
    return(rbind(mean=mean, err=err))
}

compartment.annotation = function(list.genes, cutoff.star=4)
{
    
    load(file = 'Tables_DATA/Compartment_known.Rdata')
    #known = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_knowledge_full.tsv',sep='\t',header=FALSE)
    #text = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_textmining_full.tsv',sep='\t',header=FALSE)
    #predict = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_predictions_full.tsv',sep='\t',header=FALSE)
    
    locas = c()
    for(gg in list.genes)
    {
        #cat(gg, '\n')
        if(gg=='')
        {
            locas = c(locas, NA)
        }else{
            
            gg = unlist(strsplit(as.character(gg), ';'))
            
            ll = c()
            for(gene in gg)
            {
                #gene = 'Csnk1e'
                ll = c(ll, which(known[,2]==gene & known[,7]>=cutoff.star))
                
                #ii = kk]
                #jj =
                #nn = c(nn, ii)
                #cc = c(cc, jj)
            }
            
            ii = ll[which(known[ll,4]=='Nucleus')]
            jj = ll[which(known[ll,4]=='Cytosol' | known[ll,4]=='Cytoplasm')]
            kk = ll[which(known[ll,4]=='Cytoskeleton')]
            
            if(length(ii)>0 & length(jj)==0)
            {
                locas = c(locas, 'Nucleus')
            }
            
            if(length(ii)>0 & length(jj)>0)
            {
                locas = c(locas, 'Nucleus/Cytoplasm')
            }
            
            if(length(ii)==0 & length(jj)>0 & length(kk)==0)
            {
                locas = c(locas, 'Cytoplasm')
            }
            
            if(length(ii)==0 & length(jj)>0 & length(kk)>0)
            {
                locas = c(locas, 'Cytoplasm/Cytoskeleton')
            }
            if(length(ii)==0 & length(jj)==0)
            {
                locas = c(locas, NA)
            }
            
        }
        
    }
    
    return(locas)


}

library(plotrix)
library("circular")
make_circ_coord = function(t,x,ttot=24)
{
    dt=(t[2]-t[1])*.45
    a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
    h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
    list(angles=a,heights=h)
}
circular_phase24H_histogram = function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5)
{
    #color_DHS = rgb(0.6,0,0.2)
    par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
    #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
    br=0:24
    h=hist(x, br=br,plot=FALSE)
    co=make_circ_coord(br[-1],h$counts)
    radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
}

outlier.test.ZT21 = function(x, t=c(0:15)*3, period=24)
{
  # x = data[52, ];t=c(0:15)*3;
  x = as.double(x);
  kk = which(!is.na(x)==TRUE)
  x = x[kk]
  t = t[kk]
  n=length(x)
  #mu=mean(x)
  nb.timepoints=length(x)
  #sig2=var(x)
  
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit = lm(x~c+s)
  
  #rstudent(fit)
  
  cd = cooks.distance(fit)
  cd.cutoff = 5*mean(cd)
  #cd.cutoff = 4/length(cd);
  #cd.cutoff = 4/(length(cd)-3-1)
  out.t = t[which(cd>cd.cutoff)]
  
  if(length(which(out.t==21))>0)
  {
    return(1);
  }else{
    return(0);
  }
  #a=fit$coef[2]
  #b=fit$coef[3]
  #R2=0
  #if(sig2>0) R2 =1.0-sum(fit$residuals^2)/(n-1)/sig2
}

