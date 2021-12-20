##########################################################################
##########################################################################
# Project: Microsomal proteins 
# Script purpose: analysis of microsomal proteins
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Dec 23 10:28:33 2020
##########################################################################
##########################################################################
version.DATA = 'microsomalProt'
version.analysis =  paste0(version.DATA, '_202012')
dataDir = paste0("../data/")
resDir = paste0('../results/', version.DATA)
tabDir = paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")
oldDir = '../archives/Microsomal_proteins/'

if(!dir.exists("../results/")){dir.create("../results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

########################################################
########################################################
# Section : raw data processing
# 
########################################################
########################################################
Process.rawData = FALSE
if(Process.rawData){
  start.with.Raw.Data = FALSE
  if(start.with.Raw.Data)
  {
    require(preprocessCore)
    a = read.delim('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/microsomal_proteins_new/proteinGroups.txt', 
                   header=TRUE, sep='\t',  na.strings = "Non Numérique")
    b = read.delim('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/7199-7210_microsomal_mutants/proteinGroups.txt', header=TRUE, sep='\t', 
                   na.strings = "Non Numérique")
    #mm = match(b$Majority.protein.IDs, a$Majority.protein.IDs)
    ii = grep('Ratio.H.L.normalized.ZT', colnames(a))
    jj = grep('Ratio.H.L.normalized.', colnames(b))
    aa = a[,ii]
    aa = as.matrix(log2(aa));
    aa = -aa
    bb = b[,jj]
    bb = -bb
    bb = as.matrix(bb)
    Data.control = FALSE
    if(Data.control)
    {
      boxplot(aa, ylim=c(-1, 1));abline(h=0, lwd=2.0, col='red')
      boxplot(bb, ylim=c(-10, 10));abline(h=0, lwd=2.0, col='red')
      xx = read.table('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/microsomal_proteins_new/proteinGroups_nocontaminants_log2_alldata.txt', 
                      header=TRUE, sep='\t',  na.strings = "NaN")
      yy = read.table('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/7199-7210_microsomal_mutants/Microsomal_all_values.txt', 
                      header=TRUE, sep='\t',  na.strings = "NaN")
      
      xx = as.matrix(xx[, c(1:16)]) 
      boxplot(xx, ylim=c(-10, 10));abline(h=0, lwd=2.0, col='red')
      
      yy = as.matrix(yy[, c(1:12)])
      boxplot(yy, ylim=c(-10, 10));abline(h=0, lwd=2.0, col='red')
    }
    
    bb = bb[, c(9, 12, 10, 11, 5, 8, 6, 7, 1, 4, 2, 3)]
    
    bb = data.frame(b$Gene.names, bb, stringsAsFactors = FALSE)
    colnames(bb) = c('Gene.names', 'ZT00.CryKO', 'ZT06.CryKO', 'ZT12.CryKO', 'ZT18.CryKO',
                     'ZT00.BmalWT', 'ZT06.BmalWT', 'ZT12.BmalWT', 'ZT18.BmalWT', 
                     'ZT00.BmalKO', 'ZT06.BmalKO', 'ZT12.BmalKO', 'ZT18.BmalKO')
    
    order = paste('ZT', seq(0, 46, by=3), sep='')
    cname = colnames(aa)
    cname = unlist(strsplit(as.character(cname), '[.]'))[(c(1:16)*6-1)]
    mm = match(order, cname)
    aa = aa[,mm]
    colnames(aa) = c("ZT00.WT","ZT03.WT","ZT06.WT","ZT09.WT","ZT12.WT","ZT15.WT","ZT18.WT","ZT21.WT","ZT24.WT","ZT27.WT","ZT30.WT","ZT33.WT","ZT36.WT","ZT39.WT","ZT42.WT","ZT45.WT")
    keep = c("Gene.names","Protein.IDs", "Majority.protein.IDs","Protein.names","Number.of.proteins", "Fasta.headers");
    kk = match(keep, colnames(a))
    aa = cbind(aa, a[,kk])
  }else{
    #### import processing tables (in log2)
    xx = read.table('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/microsomal_proteins_new/proteinGroups_nocontaminants_log2_alldata.txt', 
                    header=TRUE, sep='\t',  na.strings = "NaN")
    yy = read.table('/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Raw_Data/7199-7210_microsomal_mutants/Microsomal_min_one_value.txt', 
                    header=TRUE, sep='\t',  na.strings = "NaN")
    
    xx = xx[, c(1:16, 27:30)]
    yy = yy[, c(27, 9:12, 5:8, 1:4)]
    
    ### convert the H/L ratio into L/H ratio
    xx = data.frame(-as.matrix(xx[, c(1:16)]), xx[, -c(1:16)], stringsAsFactors = FALSE)
    yy = data.frame(yy[, 1], -as.matrix(yy[, c(2:13)]), stringsAsFactors = FALSE)
    colnames(xx)[1:16]  = c("ZT00.WT","ZT03.WT","ZT06.WT","ZT09.WT","ZT12.WT","ZT15.WT","ZT18.WT","ZT21.WT",
                            "ZT24.WT","ZT27.WT","ZT30.WT","ZT33.WT","ZT36.WT","ZT39.WT","ZT42.WT","ZT45.WT")
    colnames(yy) = c('Gene.names', 'ZT00.CryKO', 'ZT06.CryKO', 'ZT12.CryKO', 'ZT18.CryKO',
                     'ZT00.BmalWT', 'ZT06.BmalWT', 'ZT12.BmalWT', 'ZT18.BmalWT', 
                     'ZT00.BmalKO', 'ZT06.BmalKO', 'ZT12.BmalKO', 'ZT18.BmalKO')
    aa = xx;
    bb = yy;
    
  }
  
  find.nb.timepoints = function(data)
  {
    return(length(which(!is.na(data))));
  }
  keep = rep(NA, length=nrow(aa))
  for(n in 1:nrow(aa))
  {
    if(aa$Gene.names[n]!='')
    {
      kk = which(as.character(bb[,1])==aa$Gene.names[n])
      if(length(kk)==1) keep[n] = kk;
      if(length(kk)>1)
      {
        nn = apply(as.matrix(bb[kk, -1]), 1, find.nb.timepoints)
        kk = kk[which(nn==max(nn))]
        kk = kk[1]
        keep[n] = kk;
      } 
    }
  }
  
  xx = data.frame(aa[, c(1:16)], bb[keep, -1], aa[, -c(1:16)], stringsAsFactors = FALSE)
  #xx = xx[, -34]
  
  aa = xx;
  
  #write.table(aa,'Tables/Microsomal_proteins_WT_KO_All_L_H_ratio_log2.txt',sep='\t',quote=FALSE,col.names=TRUE, row.names = FALSE)
  
}


########################################################
########################################################
# Section : quick rhythmicity test 
# 
########################################################
########################################################
#### Start Analysis of the processed table
aa = read.table('Tables/Microsomal_Proteins_WT_KO_stat_All_L_H_ratio_log2.txt',sep='\t', header=TRUE)
source('f24_modified_1.0.r')
source('function_microsomal.R')

### Add statistics  
#aa = read.table('Microsomal_Proteins_all.txt', sep='\t', header=TRUE)
kk = which(aa$Gene.names=='')
aa$Gene.names[kk] = NA
jj = which(aa$Majority.protein.IDs=='')

source('f24_modified_1.0.r')

data = as.matrix(aa[,c(1:16)])
res = t(apply(data,1, f24_R2_alt2, t=c(0:15)*3))
res[,4] = t(apply(2^data,1, f24_R2_alt2, t=c(0:15)*3))[,4]
qv = qvals(res[,6])
res = cbind(res, qv)
#colnames(res) = paste(colnames(res), '.WT', sep='')
o = order(res[,6])
res = res[o,]
aa = aa[o,]

aa = cbind(aa, res)

source('function_microsomal.R')
res.ko = c()
for(n in 1:nrow(aa))
{
  res.ko = rbind(res.ko, model.sel.wt.ko(c(as.matrix(aa[n, c(1:16)]), as.matrix(aa[n, grep('CryKO', colnames(aa))]))))
}

colnames(res.ko) = c('prob.M1.CryKO', 'prob.M2.CryKO', 'prob.M3.CryKO', 'nb.CryKO', 'corr.CryKO')

aa = cbind(aa, res.ko)
#write.table(aa,'Tables/Microsomal_Proteins_WT_KO_stat_All_L_H_ratio_log2.txt',sep='\t',quote=FALSE, col.names=TRUE, row.names=FALSE)


CIRC = c("Per1","Per2","Per3","Cry1","Cry2","Dbp","Tef","Hlf","Nr1d1","Nr1d2","Rora","Rorb","Rorc","Arntl","Bmal1","Clock","Npas2","Bhlhe40","Bhlhe41","Cirbp","Hamp","Hamp2","Nr4a2","Tfrc","Wee1", "Por")
mm = match(CIRC, aa$Gene.names)
CIRC[which(!is.na(mm)==TRUE)]

########################################################
########################################################
# Section : Add the subcellular localization for detected proteins
# 
########################################################
########################################################
local.mouse = read.delim('uniprot-proteome_mouse_localization_GO_components.tab', sep='\t', header=TRUE)
known = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Localization_COMPARTMENTS/mouse_compartment_knowledge_full.tsv',
                   sep='\t',header=FALSE)
text = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Localization_COMPARTMENTS/mouse_compartment_textmining_full.tsv',
                  sep='\t',header=FALSE)
predict = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Localization_COMPARTMENTS/mouse_compartment_predictions_full.tsv',
                     sep='\t',header=FALSE)
load(file='/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Rdata/Uniprot_Compartments_localization_genes_all.Rdata')
#local.human = read.delim('uniprot-proteome_human_localization.tab', sep='\t', header=TRUE)
#local.rat = read.delim('uniprot-proteome_rat_localization.tab', sep='\t', header=TRUE)
secreted = read.delim('Annotations/uniprot-locations-secreted.tab', header = TRUE)
secreted = secreted$Gene.names
secreted = unique(unlist(strsplit(as.character(secreted), ' ')))

Processing.Localization.Uniprot.Compartments = FALSE
if(Processing.Localization.Uniprot.Compartments)
{
  ### list of gene names in Uniprot
  local.names = c()
  for(n in 1:nrow(local.mouse))
  #for(n in 1:100)
  {
    if(n%%1000==0) cat(n, '\n');
    gg = (local.mouse$Gene.names[n])
    gg = unlist(strsplit(as.character(gg), ' '))
    if(length(gg)==1) local.names = rbind(local.names, c(n, as.character(local.mouse$Gene.names[n]), gg));
    if(length(gg)>1)
    {
      for(g.list in gg) local.names = rbind(local.names, c(n, as.character(local.mouse$Gene.names[n]), g.list));
    }
  }
  local.names = data.frame(local.names)
  colnames(local.names) = c('index', 'gene.names.init', 'gene.names')
  
  #save(local.mouse, local.names, file='/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Rdata/Uniprot_localization_genes_names.Rdata')
  genes = unique(c(as.character(local.names[, 3]), as.character(known[,2]), as.character(text[, 2]), as.character(predict[,2])))
  
  ll = matrix(NA, nrow=length(genes), ncol=6)
  colnames(ll) = c('gene', 'uniprot.component', 'uniprot.local', 'compart.known', 'compart.text', 'compart.predict')
  
  
  #rm.curlybkt.in <- function(x) gsub("\\{.*\\}", "", x);
  rm.curlybkt.in <- function(x) 
  {
    # mm = 16;x = local.mouse[mm, 8];
    #x = as.character(x);
    x = unlist(strsplit(as.character(x), '\\. '));
    test = gsub(" \\{.*\\}", "", x);
    test =  gsub("\\{.*\\}", "", test);
    ii = grep('Note=', test)
    if(length(ii)>0) test = test[-ii];
    ii = grep('(By similarity)', test)
    if(length(ii)>0) test = test[-ii];
    #test = unlist(strsplit(as.character(test), '\\: '))
    test = gsub("\\.", "", test);
    test = gsub("SUBCELLULAR LOCATION: ", "", test);
    test = gsub("; Single-pass membrane protein", "", test);
    test = gsub("; Multi-pass membrane protein", "", test);
    removals = c("SUBCELLULAR", "LOCATION:", "", ".")
    test = test[which(is.na(match(test, removals))==TRUE)]
    return(unique(test));
  }
    
  #for(n in 1:400)
  for(n in 1:length(genes))
  {
    if(n%%1000==0) cat(n, '\n');
    
    keep = rep(NA, 6)
    gg = genes[n];
    keep[1] = gg
    
    ### extract localizaiton from Uniprot
    index = local.names[which(local.names[,3]==gg),1]
    if(length(index)>0)
    {
      subcellular = c()
      ccom = c()
      for(mm in index)
      {
        #mm = 10
        mm = as.integer(mm);
        ssub = rm.curlybkt.in(local.mouse[mm, 8])
        #print(mm);print(ssub);
        #ssub = sapply(ssub, rm.curlybkt.in, simplify = TRUE, USE.NAMES = FALSE)
        #ssub= rm.curlybkt.in(as.character(local.mouse[mm, 8]));
        #ssub = unlist(strsplit(ssub,':'))[-1]
        subcellular = c(subcellular, ssub)
        ccom = c(ccom, unlist(strsplit(as.character(local.mouse[mm, 9]), '; ')))
      }
      subcellular = unique(subcellular[which(!is.na(subcellular) & subcellular!='')])
      ccom = unique(ccom[which(ccom!='' & !is.na(ccom))])
      
      if(length(subcellular)>0){subcellular = paste(subcellular, sep='',collapse=';');keep[3] = as.character(subcellular);}
      if(length(ccom)>0){ccom = paste(ccom, sep='', collapse=';'); keep[2] = as.character(ccom);}
    }
    
    ### extract localizaiton from Compartments
    jj1 = which(known[,2]==gg)
    if(length(jj1)>0)
    {
      sub1 = paste(known[jj1, 4], sep='', collapse = ';')
      score1 = paste(known[jj1, 7], sep='', collapse = ';')
      subscore1 = paste(c(sub1, score1), sep='', collapse = '||')
      keep[4] = subscore1;
    }
    
    jj2 = which(text[,2]==gg)
    if(length(jj2)>0)
    {
      sub2 = paste(text[jj2, 4], sep='', collapse = ';')
      score2 = paste(text[jj2, 6], sep='', collapse = ';')
      subscore2 = paste(c(sub2, score2), sep='', collapse = '||')
      keep[5] = subscore2; 
    }
    
    jj3 = which(predict[,2]==gg)
    if(length(jj3)>0)
    {
      sub3 = paste(predict[jj3, 4], sep='', collapse = ';')
      score3 = paste(predict[jj3, 7], sep='', collapse = ';')
      subscore3 = paste(c(sub3, score3), sep='', collapse = '||')
      keep[6] = subscore3;
    }
    
    ll[n, ] = keep;
  }
  save(ll, file='/Users/jiwang/Proteomics_anaylysis/Microsomal_proteins/Rdata/Uniprot_Compartments_localization_genes_all.Rdata')
  
}

##########################################
### add subcellular localization for detected proteins
### Uniprot and Compartment with confidence score >=3;
##########################################
source('microsomal_functions.R')

### assign subcellular localization to detected microsomal proteins
annots = matrix(NA, nrow=nrow(aa), ncol=3)
#colnames(annots) = c('uniprot.local', 'uniprot.secreted', 'compartments.local', 'compartments.known', 'compartments.text', 'compartments.predict')
colnames(annots) = c('uniprot.secreted', 'uniprot.local',  'compartments.local')
#annots.local = c('Nucleus', 'Endoplasmic reticulum', 'Golgi apparatus', 'Cell membrane', 'Vesicle/Vacuole',
#                 'Mitochondrion', 'Extracellular space', 'Cytosol', 'Cytoskeleton', 'Peroxisome', 'Lysosome',  'Endosome')

library(plyr)
xx = count(known[,4])
xx = xx[order(-xx[,2]), ]
xx[,2] = xx[,2]/length(unique(known[,2]))

for(n in 1:nrow(aa))
#for(n in 1:200)
{
  #n = 1;
  gg = aa$Gene.names[n];
  gg = unlist(strsplit(as.character(gg), ';'))
  
  keep = rep(NA, 3);
  if(length(which(!is.na(match(gg, secreted))))>0) keep[1] = 1
  gg = gg[which(!is.na(match(gg, ll[,1])))]
  
  if(length(gg)>0)
  {
    ### localization in Uniprot with protein ID 
    index = match(gg, ll[,1])
    index = index[which(!is.na(ll[index, 2])==TRUE)];
    subs.u = paste(ll[index, 3], sep='', collapse = ';');
    keep[2] = find.uniprot.localization(subs.u)
    
    ### localization from Compartment with gene names
    jj1 = c();
    for(g.test in gg) jj1 = c(jj1, which(known[,2]==g.test));
    jj1 = unique(jj1);
    if(length(jj1)>0)
    {
      jj = jj1[which(known[jj1, 7]>=3)] 
      subs.c = unique(known[jj, 4])
      subs.c = paste(subs.c, sep='', collapse = ';');
      keep[3] =  find.uniprot.localization(subs.c)
    }
  }
  annots[n, ] = keep;
}

aa = data.frame(aa, annots, stringsAsFactors=FALSE)

length(which(!is.na(aa$uniprot.local)))
length(which(!is.na(aa$compartments.local)))
length(which(!is.na(aa$uniprot.local) |!is.na(aa$compartments.local)))

#colnames(aa)[30:31] = c('subcellular.localization', 'GO.cell.components')
write.table(aa,'Tables/Microsomal_Proteins_L_H_ratio_log2__all_stat_localization_Uniprot_COMPARTMENTS.txt',sep='\t',quote=FALSE,col.names=TRUE)

###########
#### Plot gene profiles
###########
#### compare profiles of wt and cryko nuclear proteins
pdf('Plots/Microsomal_profiles_all_genes_not_centered_WT_KO_log2_with_localization.pdf',width=15,height=8)
#for(n in 1:20)
for(n in 1:nrow(aa))
{
	if(aa$nb.timepoints[n]>=8) ## select proteins with at least 4 time points
	{
		name = aa$Gene.names[n]
		y0 = as.numeric(aa[n,c(1:16)])
		#y1 = as.numeric(aa[n,c(1:16)])
		#y0 = (y0-mean(y0[which(!is.na(y0)==TRUE)]))
				
		time = c(0:15)*3
		lims = range(as.numeric(aa[n, 1:28]), na.rm = TRUE)
		
		plot(time, y0, ylim=lims, col='blue', type='b',lwd=2.0, ylab='Protein abundance (log2)', xlab='ZT[h]',lty=1, 
		     main=paste(name, " phase=", signif(aa$phase[n],d=2), ", pval=", signif(aa$pval[n],d=2), ", FDR=", signif(aa$qv[n],d=2), 
		                '\n', 'Uniprot: ', aa$uniprot.local[n], '  Secreted: ', aa$uniprot.secreted[n], '\nCOMPART: ', aa$compartments.local[n],  sep=''))
		points(c(0:3)*6, as.numeric(aa[n, c(17:20)]), type='b', lwd=2.0, col='red')
		points(c(0:3)*6, as.numeric(aa[n, c(21:24)]), type='b', lwd=2.0, col='green')
		points(c(0:3)*6, as.numeric(aa[n, c(25:28)]), type='b', lwd=2.0, col='orange')
		abline(h=0.0,col='darkgray',lwd=3.0)
		legend('topright', legend = c('WT', 'Cry DKO', 'Bmal1 WT', 'Bmal1 KO'), col='black', fill = c('blue', 'red', 'green', 'orange'), bty='n');
		
	}
}

dev.off()


########################################################
########################################################
# Section : Quality control and further rhythmicity analysis
# 
# reanalysis were done after the location annotation
########################################################
########################################################
microsomal = read.table(file = paste0(dataDir, 'microProt_L.H.ratio_log2_all_stat_localization.Uniprot.COMPARTMENTS.txt'),
                        sep='\t', header = TRUE)
data = as.matrix(microsomal[, c(1:28)])

Check.QCs.with.plots = FALSE
if(Check.QCs.with.plots){
  pdf.name = paste0(resDir, "/Boxplots_all_samples_WT_KO.pdf")
  pdf(pdf.name, width=15, height=10)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(6,3,2,0.8)+0.1, tcl = -0.3)
  
  boxplot(data, las = 2)
  
  dev.off()
  
  #pairs(data, na.action = stats::na.pass)
  
  pdf.name = paste0(resDir, "/Correlation_replicates_WT_KO.pdf")
  pdf(pdf.name, width=8, height=6)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
  
  par(mfrow=c(3, 4))
  cex = 0.02;
  for(n in 1:8)
  {
    #print(c((n-1)*3, (n+7)*3))
    plot(data[, c(n, (n+8))], cex=cex, xlab=paste('ZT', (n-1)*3, '-Rep1', sep = ''),
         ylab=paste('ZT', (n-1)*3,  '-Rep2', sep=''), xlim=c(-4, 4), ylim=c(-4, 4), main=NA);
    jj = which(!is.na(data[,n]==TRUE & !is.na(data[, (n+8)])))
    test = (cor.test(data[jj,n], data[jj, (n+8)], alternative='greater'))
    
    text(-2, 3, paste('R = ', signif(cor(data[,n], data[, (n+8)], use="na.or.complete"), d=2), sep=''), col='red')
    
    abline(0, 1, lwd=1.5, col='red')
  }
  
  for(n in 1:4)
  {
    plot(data[, c((n*2-1), (n+16))], cex=cex, xlab=paste('WT-ZT', (n-1)*6, '-Rep1', sep = ''),
         ylab=paste('CryDKO-ZT', (n-1)*6, sep=''), xlim=c(-4, 4), ylim=c(-4, 4), main=NA);
    #jj = which(!is.na(data[,n]==TRUE & !is.na(data[, (n+8)])))
    #test = (cor.test(data[jj,n], data[jj, (n+16)], alternative='greater'))
    text(-2, 3, paste('R = ', signif(cor(data[,(n*2-1)], data[, (n+16)], use="na.or.complete"), d=2), sep=''), col='red')
    abline(0, 1, lwd=1.5, col='red')
  }
  #cd = cooks.distance(fit)
  dev.off()
  
}

##########################################
# Rhythmicity analysis
##########################################
source('f24_modified_1.0.r')

data = as.matrix(microsomal[,c(1:16)])

period = 12
res = t(apply(data,1, f24_R2_alt2, t=c(0:15)*3, period=period))
res[,4] = t(apply(2^data,1, f24_R2_alt2, t=c(0:15)*3, period=period))[,4]
qv = qvals(res[,6])
res = cbind(res, qv)

colnames(res) = paste0(colnames(res), '.WT.12hRhythmicity')

jj = match(c("nb.timepoints", "mean", "amp", "relamp", "phase", "pval", "qv"), colnames(microsomal))
stats = microsomal[, jj]
xx = microsomal[, -jj]
colnames(stats) = paste0(colnames(stats), '.WT.24hRhythmicity')

mcm = data.frame(xx, stats, res, stringsAsFactors = FALSE)

write.csv(mcm, file = paste0(tabDir, 'microsomalProt_processed.log2_localizationAnnot_rhythmicity.csv'), row.names = FALSE)
saveRDS(mcm, file = paste0(RdataDir, 'microsomalProt_processed.log2_localizationAnnot_rhythmicity.rds'))

Make.plots.for.rhythmicity = FALSE
if(Make.plots.for.rhythmicity){
  qq = c(0:100)/100
  nb.rhythmic = c()
  for(n in 1:length(qq))
  {
    cutoff.qq = qq[n]
    nb.rhythmic = c(nb.rhythmic,length(which(microsomal$qv<cutoff.qq)))
  }
  
  pdf('Plots/Nb_rhythmic_Proteins_FDR.pdf', width=2., height=2.)
  par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  
  plot(qq, nb.rhythmic+0.01, type='l', lty=1, lwd=1.2, main=NA, log='y',ylim=c(1,5000), xlab=NA, ylab=NA, col=c('black'), axes=FALSE)
  #points(qq, nb.rhythmic.robles+0.01, type='l', lty=1, lwd=2, col='green')
  abline(v=0.15,col='darkblue',lwd=1.5,lty=2)
  abline(v=0.05,col='darkblue',lwd=1.5,lty=2)
  abline(v=0.1,col='darkblue',lwd=1.5,lty=2)
  axis(1,at=seq(0, 1.0, by=0.2),cex.axis =1.0)
  axis(2, at= c(1, 5, 50, 500, 5000), las=1,cex.axis = 1.0)
  box()
  #abline(v=0.2,col='darkgreen',lwd=2,lty=2)
  #abline(v=0.1,col='darkgreen',lwd=2,lty=2)
  #legend('topright', legend = c('Mauvoisin','Robles'), lty=c(1,1), cex=0.7,col = c('blue', 'green'), border = NA, bty = 'n')
  dev.off()
  
  source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
  
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
    #color_DHS = rgb(0.6,0,0.2);x=phases;
    par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
    #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
    br=0:24
    h=hist(x, br=br,plot=FALSE)
    co=make_circ_coord(br[-1],h$counts)
    radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2, main=NA, rp.type='p', poly.col=color_hist)
  }
  
  #### All detected rhythmic proteins
  jj = which(microsomal$qv<0.15)
  phases = as.numeric(microsomal$phase[jj])
  amps = microsomal$amp[jj]
  
  pdf('Plots/Phase_distribution_all_rhythmic_detected.pdf', width=2.2, height=2.2)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  
  #circular_phase24H_histogram(phases, col=rgb(0.6,0,0.2), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
  hist(phases, breaks=c(0:12)*2, col='gray')
  
  dev.off()
  
  pdf('Plots/Amplitudes_distribution_all_rhythmic_detected.pdf', width=2.2, height=2.2)
  par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
  breaks=c(0:10)/2;
  h = hist(amps, breaks=breaks, plot=FALSE)
  y = h$counts
  y[which(y<=0)] = 0.5;
  lwd = 1.2
  plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd)
  axis(1,at=seq(0, 5, by=1),cex.axis =1.0)
  axis(2, at= c(1, 2, 5, 10, 20, 50, 100, 200), las=1,cex.axis = 1.0)
  lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd)
  lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd)
  lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd)
  lines(h$breaks, rep(0,length(h$breaks)), type='S')
  invisible(h)
  
  dev.off()
  
  #### Check phases and amplitudes
  kk = which(microsomal$qv<0.15)
  genes = microsomal$Gene.names[kk]
  genes = unlist(strsplit(as.character(genes), ';'))
  genes = genes[which(!is.na(genes)==TRUE)]
  write.table(genes, file='Tables/rhythmic_microsomal_proteins_4functional_analysis.txt', sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
  
}

########################################################
########################################################
# Section : Combine microsomal proteins with mRNA, total proteins, nuclear proteins.
# 
########################################################
########################################################
microsomal = read.table('Tables/Microsomal_Proteins_L_H_ratio_log2__all_stat_localization_Uniprot_COMPARTMENTS.txt',sep='\t', header = TRUE);

microsomal.names = c()
for(n in 1:nrow(microsomal))
{
  if(!is.na(microsomal$Gene.names[n])==TRUE)
  {
    gg = unlist(strsplit(as.character(microsomal$Gene.names[n]), ';'));
    for(g.test in gg) microsomal.names = rbind(microsomal.names, c(n, as.character(microsomal$Gene.names[n]), g.test));
  }
}
microsomal.names = data.frame(microsomal.names, stringsAsFactors = FALSE);
colnames(microsomal.names) = c('index', 'Gene.names', 'gene');

#### Import mRNAs, total and nuclear proteins
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/mRNA_WT_RF_total_RNA_seq_Cedric.Rdata')
mrna.names = data.frame(c(1:nrow(mrna)), mrna$gene, mrna$nb.timepoints.mrna, mrna$pval.mrna, stringsAsFactors=FALSE)

total = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Table_total_proteins_PNAS.txt',sep='\t', header=TRUE)
total.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Table_total_proteins_index_names.txt', 
                         sep='\t', header=TRUE, as.is=c(2))
total.names = data.frame(total.names, total[total.names[,1], c(17, 23)], stringsAsFactors=FALSE)

nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', 
                     sep='\t', header=TRUE, as.is=c(17:20))
nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',
                           as.is=c(2,3), header=TRUE, sep='\t')
nuclear.names = data.frame(nuclear.names[, c(1,3)], nuclear$nb.timepoints[nuclear.names[,1]], nuclear$pval[nuclear.names[,1]], stringsAsFactors = FALSE)


index.mrna = rep(NA, nrow(microsomal))
index.total = rep(NA, nrow(microsomal));
index.nuclear = rep(NA, nrow(microsomal));

for(n in 1:nrow(microsomal))
{
  if(!is.na(microsomal$Gene.names[n])==TRUE)
  {
    ii = which(microsomal.names[,1]==n)
    ggene = microsomal.names[ii,3]
    
    mapping.total = c()
    mapping.mrna = c()
    mapping.nuclear = c()
    for(gene in ggene)
    {
      mapping.total = c(mapping.total, which(total.names[,2]==gene))
      mapping.mrna = c(mapping.mrna, which(mrna.names[,2]==gene))
      mapping.nuclear = c(mapping.nuclear, which(nuclear.names[,2]==gene))
    }
    mapping.total = unique(mapping.total)
    mapping.mrna = unique(mapping.mrna)
    mapping.nuclear = unique(mapping.nuclear)
    
    ## mapping total
    mapping = mapping.total
    if(length(mapping)>=1)
    {
      mapping = mapping[which(total.names[mapping,4]==min(total.names[mapping,4]))];
      index.total[n] = total.names[mapping[1], 1];
    }
    ## mapping mrna
    mapping = mapping.mrna
    if(length(mapping)>=1)
    {
      mapping = mapping[which(mrna.names[mapping,4]==min(mrna.names[mapping,4]))];
      index.mrna[n] = mrna.names[mapping[1], 1];
    }
    ## mapping pol2
    mapping = mapping.nuclear
    if(length(mapping)>=1)
    {
      mapping = mapping[which(nuclear.names[mapping,4]==min(nuclear.names[mapping,4]))];
      index.nuclear[n] = nuclear.names[mapping[1], 1];
    }
    
  }
}

mapping.total = total[index.total,]
mrna$qv.mrna = qvals(mrna$pval.mrna)
mapping.mrna = mrna[index.mrna,]
mapping.nuclear = nuclear[index.nuclear,]

colnames(mapping.mrna)[1] = 'gene.mrna'

mapping.total = mapping.total[, c(29, 1:16, 18:24)]
colnames(mapping.total)[c(1,18, 24)] = c('gene.total', 'nb.timepoints.total', 'qv.total')

#mapping.total = mapping.total[, c(1:30)]
#colnames(mapping.total)[29] = 'gene.total'
mapping.nuclear = mapping.nuclear[, c(20,1:16, 22:28)]
colnames(mapping.nuclear)[1] = 'gene'
colnames(mapping.nuclear) = paste(colnames(mapping.nuclear), '.nuclear', sep='')
#colnames(mapping.pol2)[1] = 'gene.pol2'
#colnames(mapping.mrna)[c(32, 38)] = c('qv.mrna', 'qv.12h.mrna')
#xx = mapping.total[, c(1:16,17,19:24,29,31:35)]
#mapping.total = xx
#colnames(mapping.total)[c(17,23,29)] = c('nb.timepoints.total', 'qv.total', 'qv.12h.total')
microsomal.all = data.frame(microsomal, mapping.mrna, mapping.total, mapping.nuclear, stringsAsFactors=FALSE)

HPA.predicted.secretome = FALSE
if(HPA.predicted.secretome)
{
  hpa = read.delim('Annotations/Secreted_proteins_predicted_HPA.tab', header = TRUE)
  xx = sapply(hpa[,1], first.upper)
  hpa[, 1] = xx
  hpa.secreted = rep(NA, nrow(microsomal))
  
  mm = match(microsomal.names[,3], hpa[,1])
  kk = which(!is.na(mm)==TRUE)
  index = as.integer(unique(microsomal.names[kk, 1]))
  hpa.secreted[index] = 1
  microsomal = data.frame(microsomal[, c(1:44)], hpa.secreted,  microsomal[, c(45:47)], stringsAsFactors = FALSE)
  microsomal.all = data.frame(microsomal.all[, c(1:44)], hpa.secreted,  microsomal.all[, c(45:115)], stringsAsFactors = FALSE)
}

save(microsomal, microsomal.names, microsomal.all, file='Rdata/Table_microsomal_log2_L_H_names_mRNA_total_nuclear.Rdata')
#write.table(microsomal.all, 'Tables/Table_microsomal_mRNA_total_nuclear.txt', sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)

#### Check some numbers
nrow(microsomal.all)

length(which(microsomal.all$qv<0.15))
length(which(microsomal.all$qv.mrna<0.1))/length(which(!is.na(microsomal.all$qv.mrna)))
length(which(microsomal.all$qv.total<0.25))/length(which(microsomal.all$qv.total<=1))
length(which(microsomal.all$pval.nuclear<0.05))/length(which(!is.na(microsomal.all$pval.nuclear)))

### overlap between microsomal and mRNAs
length(which(microsomal.all$qv.mrna<0.1))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.mrna<=1))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.mrna<0.1))
length(which(microsomal.all$qv<0.15 & microsomal.all$uniprot.secreted==1))
length(which(microsomal.all$qv<0.15 & microsomal.all$uniprot.secreted==1 & microsomal.all$qv.mrna<0.1))
length(which(microsomal.all$qv<0.15 & microsomal.all$hpa.secreted==1))
length(which(microsomal.all$qv<0.15 & microsomal.all$hpa.secreted==1 & microsomal.all$qv.mrna<0.1))

### overlap between microsomal and total proteins
length(which(microsomal.all$qv.total<0.25))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.total<=1))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.total<=0.25))

### overlap between microsomal and nuclear 
length(which(microsomal.all$qv.nuclear<0.05))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.nuclear<=1))
length(which(microsomal.all$qv<0.15 & microsomal.all$qv.nuclear<=0.05))

#### compare profiles of protein profiles in microsomal, total, nuclear and mRNA
pdf('Plots/Microsomal_profiles_all_genes_centered_WT_KO_log2_with_localization_mRNA_total_nuclear.pdf',width=15,height=8)
#for(n in 1:20)
aa = microsomal.all
for(n in 1:nrow(aa))
{
  if(aa$nb.timepoints[n]>=8) ## select proteins with at least 8 time points
  {
    name = aa$Gene.names[n]
    y0 = as.numeric(aa[n,c(1:16)])
    y1 = as.numeric(aa[n, grep('mRNA.ZT', colnames(aa))])
    y2 = as.numeric(aa[n, grep('total', colnames(aa))[2:17]])
    y3 = as.numeric(aa[n, grep('nuclear', colnames(aa))[2:17]])
    for(jj in c(0:3)) 
    {
      eval(parse(text=paste('y = y', jj, sep = '')));
      y = y - mean(y, na.rm = TRUE);
      eval(parse(text=paste('y', jj, ' = y;', sep = '')));
    }
      
    #y1 = as.numeric(aa[n,c(1:16)])
    #y0 = (y0-mean(y0[which(!is.na(y0)==TRUE)]))
    
    time = c(0:15)*3
    lims = range(c(y0, y1, y2, y3), na.rm = TRUE)
    
    plot(time, y0, ylim=lims, col='blue', type='b',lwd=2.0, ylab='Protein abundance (log2)', xlab='ZT[h]',lty=1, 
         main=paste(name, " phase=", signif(aa$phase[n],d=2), ", pval=", signif(aa$pval[n],d=2), ", FDR=", signif(aa$qv[n],d=2), 
                    '\n', 'Uniprot: ', aa$uniprot.local[n], '  Secreted: ', aa$hpa.secreted[n], '\nCOMPART: ', aa$compartments.local[n],  sep=''))
    points(c(0:11)*4, y1, type='b', lwd=2.0, col='green')
    points(time, y2, type='b', lwd=2.0, col='orange', lty=2)
    points(time, y3, type='b', lwd=2.0, col='darkgray', lty=2)
    #abline(h=0.0,col='darkgray',lwd=3.0)
    legend('topright', legend = c('Microsomal', 'mRNA', 'total', 'nuclear'), col='black', fill = c('blue', 'green', 'orange', 'darkgray'), bty='n');
    
  }
}

dev.off()

########################################################
########################################################
# Section : Protein complexes analysis
# 
########################################################
########################################################
Analysis.protein.complexes = FALSE
if(Analysis.protein.complexes)
{
  source('proteomics_functions.R')
  source('microsomal_functions.R')
  microsomal = readRDS(file = paste0(RdataDir, 'microsomalProt_processed.log2_localizationAnnot_rhythmicity.rds'))
  microsomal$nb.timepoints = microsomal$nb.timepoints.WT.24hRhythmicity
  
  # prepare the detected genes in microsomal
  #microsomal = mcm
  
  microsomal.names = c()
  for(n in 1:nrow(microsomal))
  {
    if(!is.na(microsomal$Gene.names[n])==TRUE)
    {
      gg = unlist(strsplit(as.character(microsomal$Gene.names[n]), ';'));
      for(g.test in gg) microsomal.names = rbind(microsomal.names, c(n, as.character(microsomal$Gene.names[n]), g.test));
    }
  }
  microsomal.names = data.frame(microsomal.names, stringsAsFactors = FALSE);
  colnames(microsomal.names) = c('index', 'Gene.names', 'gene');
  
  # prepare the CORUM protein complexes annotation
  # annot.corum = Processing.CORUM.all.Complexes.add.manually()
  
  ##########################################
  # check detected subunits for protein complexes
  ##########################################
  Corum.table.generate = FALSE
  if(Corum.table.generate)
  {
    load(file= paste0(oldDir, 'Rdata/Annotation_Corum_all_with_Manual.Rdata'))
  
    ii = which(annot.corum$Synonyms=='')
    annot.corum$Synonyms[ii] = NA
    colnames(annot.corum)[5] = 'subunits.UniProt.IDs'
    
    microsomal.names = data.frame(as.integer(microsomal.names[,1]), microsomal.names[, c(2,3)], stringsAsFactors = FALSE)
    
    options(warn=1)
    detected = c()
    index = c()
    nb.detected = c()
    for(n in 1:nrow(annot.corum))
    {
      cat(n, '\n')
      test = annot.corum$subunits[n]
      test = unlist(strsplit(as.character(test),","))
      
      gene.detected = c()
      index.detected = c()
      
      for(ttest in test)
      {
        kk = which(microsomal.names[,3]==ttest)
        if(length(kk)>=1) 
        {
          #print(ttest)
          iindex = as.integer(microsomal.names[kk, 1])
          iindex = iindex[which(microsomal$nb.timepoints[iindex]==max(microsomal$nb.timepoints[iindex]))][1]
          index.detected = c(index.detected, iindex)
          gene.detected = c(gene.detected, ttest)
        }
        
      }
      index.detected = unique(index.detected)
      nb.detected = c(nb.detected, length(index.detected))
      if(length(index.detected)>=1)
      {
        detected = c(detected, paste(gene.detected, sep='', collapse=','))
        index = c(index, paste(index.detected, collapse=','))
      }else{
        detected = c(detected, NA)
        index = c(index, NA)
      }
    }
    
    annot.corum = data.frame(annot.corum, detected, index, nb.detected, stringsAsFactors=FALSE)
    annot.corum$percent.detected = as.numeric(annot.corum$nb.detected)/as.numeric(annot.corum$nb.subunits)
    colnames(annot.corum)[8:11] = c('subunits.detected', 'index.detected', 'nb.detected', 'percent.detected')
    #kk = which(annot.corum$nb.detected>1)
    #annot.corum = annot.corum[kk,]
    save(annot.corum, file=paste0(RdataDir, 'Annotation_Corum_all_with_Manual_detected_microsomal.Rdata'))
    
  }
  
  ##########################################
  # test rhythmicity and similarity of detected protein complexes
  ##########################################
  load(file= paste0(RdataDir, 'Annotation_Corum_all_with_Manual_detected_microsomal.Rdata'))
  source('f24_modified_1.0.r')
  source('proteomics_functions.R')
  
  SVD.Analysis.protein.complexes = FALSE
  if(SVD.Analysis.protein.complexes)
  {
    annot = annot.corum
    kk = which(annot$nb.subunits>0)
    annot = annot[kk,]
    
    kk = which(annot$nb.detected>1)
    annot = annot[kk,]
        
    pdfname = paste0(resDir, '/SVD_PC_plot_microsomal.pdf')
    res = matrix(NA, ncol=10, nrow=nrow(annot))
    
    
    ### main function of protein complex svd test
    #res = statistics.complexes.svd(annot[ii,], res[ii,], pdfname)
    res = statistics.complexes.svd(annot=annot, prot.data = microsomal, res = res, pdfname = pdfname, TEST = FALSE)
    
    res = data.frame(annot, res, stringsAsFactors=FALSE)
    
    kk = which(res$nb.subunits.quantified>=2)
    res = res[kk,]
    dim(res)
    
    save(annot, res, file= paste0(RdataDir, 'Annotation_Corum_all_with_Manual_detected_microsomal_SVD.Rdata'))
    
    ### reduce the redundancy of PC table
    load(file= paste0(RdataDir, 'Annotation_Corum_all_with_Manual_detected_microsomal_SVD.Rdata'))
    res = reduce.redundancy.protein.complex(res)
    
    ### calculate svd statistics using precomputed background
    load(file='../archives/Rdata_nuclearProt_epfl/Background_PC_SVD_v3.Rdata') 
    nb.samples = 1000
    bg0 = rep(NA, nrow(res))
    bg1 = rep(NA, nrow(res))
    bg2 = rep(NA, nrow(res))
    pval.svd = rep(NA, nrow(res))
    bg.cutoff = rep(NA, nrow(res))
    
    for(n in 1:nrow(res))
    {
      kk = which(bg.samples==res$nb.subunits.quantified[n])
      test = bg.samples[kk, c(2:(nb.samples+1))]
      
      bg0[n] = mean(test)
      bg1[n] = mean(test) + sd(test)
      bg2[n] = mean(test) + 2.33*sd(test)
      pval.svd[n] = length(which(test>=res$svd.1st.component.p1[n]))/nb.samples
      
      ttest = test[order(-test)]
      bg.cutoff[n] = ttest[50]
      
    }
    
    res$bg0 = bg0
    res$bg1 = bg1
    res$bg2 = bg2
    res$pval.svd = pval.svd
    res$bg.cutoff = bg.cutoff
    res$qv.svd = qvals(res$pval.p1)
    
    length(which(res$pval.svd<0.05))
    length(which(res$pval.svd<0.05 & res$qv.svd < 0.15))
    
    source('proteomics_functions.R')
    #kk = which(res$pval.svd < 0.05 & res$pval.p1 < 0.05)
    kk = c(1:nrow(res))
    #source('function_microsomal.R')
    #source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
    Plots.complexes.subunits(annot = res, nn = kk, prot.data = microsomal, 
                             pdfname = paste0(resDir, '/PCs_Corum_subunits_Selected_SVD_pval_0.05.pdf'))
    
    write.csv(res, file = paste0(tabDir, 'microsomalProt_protComplex_SVDanalysis_v2.csv'), 
              row.names = FALSE)
    
  }
  
}
