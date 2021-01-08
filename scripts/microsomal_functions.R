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

first.upper = function(character)
{
  test = character
  test = tolower(test)
  test = c(toupper(unlist(strsplit(as.character(test),""))[1]), unlist(strsplit(as.character(test),""))[-1])
  test = paste(test, collapse = '')
  return(test)
}

stadardization.nona = function(x)
{
  xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
  return(xx)
}

##########################################
# protein complex rhythmicity analysis 
##########################################
statistics.complexes.svd = function(annot, res, pdfname='myplots/SVD_PC_plot.pdf')
{	
  pdfname = pdfname;
  colnames(res) = c('nb.subunits.quantified', 'd.entropy', 'p1','p2', 'amp.p1', 'phase.p1','pval.p1', 'amp.p2', 'phase.p2','pval.p2')
  load(file='Rdata/Table_microsomal_log2_L_H_names_mRNA_total_nuclear.Rdata')
  nuclear = microsomal;
  ## ii = grep('RNA polymerase II holoenzyme complex', annot[,2])[1];vect = annot[ii,]
  ##
  #annot = annot[c(1:10), ]
  #pdfname = pdfname;
  #print(vect[2],'\n')
  #cutoff.rhy = 0.01
  #cutoff.sim = 0.80
  cutoff.nb.timepoints = 8
  pdf(pdfname, width=14, height=12)
  for(n in 1:nrow(annot))
  {
    cat(n, '\n');
    vect = annot[n,]
    index = vect[which(names(vect)=='index.detected')];
    subunits = vect[which(names(vect)=='subunits.detected')];
    index = as.numeric(unlist(strsplit(as.character(index), ',')))
    subunits = unlist(strsplit(as.character(subunits), ','))
    subunits = subunits[which(nuclear$nb.timepoints[index]>=cutoff.nb.timepoints)]
    index = index[which(nuclear$nb.timepoints[index]>=cutoff.nb.timepoints)]
    
    res[n,1] = length(index)
    if(length(index)>=2)	
    {
      test = matrix(data = NA, nrow = length(index), ncol = 16)
      for(kk in 1:length(index))
        test[kk, ] = standadization.nona.impute(as.numeric(nuclear[index[kk], c(1:16)]))
      rownames(test) = subunits
      ###SVD
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
      plot(c(1:length(pl)), pl, type='b', pch=16, lwd=2.0,cex=1.0, ylab='% variance explained by each component', xlab='Index of components', ylim=c(0,1),main=paste('d= ', signif(dd, d=2),sep=''))
      abline(v=(which(pl>0.7/length(pl))[length(which(pl>0.7/length(pl)))]+0.2), col='red',lwd=2.5)
      plot(c(0:15)*3, apply(ss$u, 2, mean)[1]*ss$v[,1], type='b',lwd=2.0,col='darkblue',ylab='1st component',main=paste('amp=', signif(stat1[3],d=2), ', phase= ', signif(stat1[5], d=2), ', pval=', signif(stat1[6], d=3)))
      abline(h=0, col='gray',lwd=2.0)
      plot(c(0:15)*3, apply(ss$u, 2, mean)[2]*ss$v[,2], type='b',lwd=2.0, col='darkgreen', ylab='2nd component', main=paste('amp=', signif(stat2[3],d=2), ', phase= ', signif(stat2[5], d=2), ', pval=', signif(stat2[6], d=3)))
      abline(h=0, col='gray',lwd=2.0)
    }
  }
  
  dev.off()
  
  return(res)
}


Plots.complexes.subunits = function(annot, nn, pdfname)
{
  #table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v2.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
  load(file='Rdata/Table_microsomal_log2_L_H_names_mRNA_total_nuclear.Rdata')
  nuclear = microsomal;
  
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
        text.col = 'black';
      
        type.plot='b';
        points(c(0:15)*3, (test[i,]), type=type.plot, cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1)
        #points(c(0:15)*3, test[i,], type='l', cex=2.0, lwd=2.0, col = rainbow[(i.rel-1)%%nrow(test)+1], lty = 2)
        legend(x = 0.75*xlim[2],y = ylim[2]*(1-i.rel*min(1.0/20,(1.0/nrow(test)))), legend = paste(genes[i],', pval= ', pvals[i], ', amp= ', amps[i],sep=''), text.col=text.col, col =  rainbow[(i.rel-1)%%nrow(test)+1], lty = (i.rel-1)%/%nrow(test)+1, bty = 'n' )
      }
    }
  }
  dev.off()
  
}

