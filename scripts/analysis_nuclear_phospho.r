####################################################################################
########################################## Pre-processing the data
####################################################################################
require(preprocessCore)
aa = read.csv('DATAPHOSPHONUCLEI.csv',header=TRUE)  ###### all proteins measured with MS
ii = grep('Ratio.L.H.normalized', colnames(aa))

data = aa[,ii]
#data = data
data = data[,c(13:28, 9:12, 5:8, 1:4)]
colnames(data)[1:28] = c("ZT00.WT","ZT03.WT","ZT06.WT","ZT09.WT","ZT12.WT","ZT15.WT","ZT18.WT","ZT21.WT","ZT24.WT","ZT27.WT","ZT30.WT","ZT33.WT","ZT36.WT","ZT39.WT","ZT42.WT","ZT45.WT","ZT00.Cry.KO","ZT06.Cry.KO","ZT12.Cry.KO","ZT18.Cry.KO","ZT00.Bmal.WT","ZT06.Bmal.WT","ZT12.Bmal.WT","ZT18.Bmal.WT","ZT00.Bmal.KO","ZT06.Bmal.KO","ZT12.Bmal.KO","ZT18.Bmal.KO")
kk = match(c('Intensity','Modified.sequence', 'Protein','Position' ),colnames(aa))
bb = aa[, c(1:6,8:10, kk)]

nuclear.phos = data.frame(data, bb, stringsAsFactors=FALSE)

gene.name=nuclear.phos$Gene.names
keep = c()
trace = c()
for(n in 1:length(gene.name))
{
	trace = c(trace, as.character(gene.name[n]))
	test =  unlist(strsplit(as.character(gene.name[n]), ";"))
	test = unique(test)
	if(length(test)<1) { print("NO name"); keep = c(keep, "");}
	if(length(test)==1) 
	{
		keep = c(keep, test)
	}
	if(length(test)>1)
	{
		#print(n)
		test = paste(test, sep=" ",collapse = ";") 
		keep = c(keep, test)
	}
}
nuclear.phos$gene = keep
phos.n = nuclear.phos
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
res = t(apply(as.matrix(phos.n[,c(1:16)]), 1, f24_R2_alt2, t=c(0:15)*3))
qv = qvals(res[,6])
res = cbind(res, qv)
phos.n = data.frame(phos.n, res, stringsAsFactors=FALSE)
o1 = order(phos.n$pval)
phos.n = phos.n[o1,]

index = c(1:nrow(phos.n))
gene = phos.n$gene

iidex = c()
gg = c()
for(n in 1:length(gene))
{
	test = gene[n]
	ttest = unlist(strsplit(as.character(test), ";"))
	if(length(ttest)>0)
	{
		gg = c(gg, ttest)
		iidex = c(iidex, rep(index[n], length(ttest)))
	}
}
gg.u = unique(gg)
index = c()
for(n in 1:length(gg.u))
{
	kk = which(gg==gg.u[n])
	index = c(index, paste(iidex[kk], sep='',collapse=';'))
}

phos.n.names = data.frame(index, gg.u,stringsAsFactors=FALSE)
colnames(phos.n.names) = c('index', 'gene.name.unique')

mm = match(phos.n.names[,2], nuclear.names[,3])
length(which(is.na(mm)))
phos.n.names[which(is.na(mm)==TRUE), ]

save(phos.n, phos.n.names, file='Rdata/Phos_nuclear_all_data.Rdata')

#####
##### Add the nuclear proteins
#####

nuclear = read.table('Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')

mapping = c()
mapping.gene = c()
miss = c()
for(n in 1:nrow(phos.n))
{
	gg = phos.n$gene[n]
	gg = unlist(strsplit(as.character(gg), ";"))
	if(length(gg)>0)
	{
		kk = c()
		for(ggg in gg) kk = c(kk, nuclear.names[which(nuclear.names[,3]==ggg), 1])
		kk = unique(kk)
		if(length(kk)==0) 
		{
			print('No Mapping');
			miss = c(miss, n)
			mapping = rbind(mapping, rep(NA, 16));
			mapping.gene = c(mapping.gene, NA)
		}
		if(length(kk)==1) 
		{
			mapping = rbind(mapping, nuclear[kk, c(1:16)]);
			mapping.gene = c(mapping.gene, nuclear$Gene.names[kk])
		}
		if(length(kk)>1)
		{
			kk = kk[which(nuclear$nb.timepoints[kk]==max(nuclear$nb.timepoints[kk]))][1];
			mapping = rbind(mapping, nuclear[kk, c(1:16)]);
			mapping.gene = c(mapping.gene, nuclear$Gene.names[kk])
		}
	}else{
		mapping.gene = c(mapping.gene, NA)
		mapping = rbind(mapping, rep(NA, 16))
	}
}
colnames(mapping) = paste('ZT', c(0:15)*3, '.nuclear', sep='')

index.miss.mapping = miss

phos.nuclear = data.frame(phos.n, mapping.gene, mapping)

save(phos.n, phos.n.names, phos.nuclear, index.miss.mapping, file='Rdata/Phos_Nuclear_Mapping_all.Rdata')

#######
####### plot the phospho and nuclear proteins:
#######
load(file='../Rdata/Phos_Nuclear_Mapping_all.Rdata')
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')

pdf('myplots/Phosphos_Nuclear_profiles_centered.pdf',width=10,height=8)
for(n in 1:nrow(phos.n.names))
{
	index = phos.n.names[n,1]
	gene = phos.n.names[n, 2]
	iindex = as.numeric(unlist(strsplit(as.character(index), ";")))
	iindex = iindex[which(phos.nuclear$nb.timepoints[iindex]>=6)]
	if(length(iindex)>0)
	{
		for(index in iindex)
		{
			y1 = as.matrix(phos.nuclear[index,c(1:16)])
			pval = paste(signif(phos.nuclear$pval[index], d=2), sep='',collapse=';')
			phase = paste(signif(phos.nuclear$phase[index], d=2), sep='',collapse=';')
			y0 = as.numeric(phos.nuclear[index[1], c(51:66)])
			y2 = as.matrix(phos.nuclear[index,c(17:20)])
			y3 = as.matrix(phos.nuclear[index,c(21:24)])
			y4 = as.matrix(phos.nuclear[index,c(25:28)])
			#y0 = t(apply(y0, 1, centering.nona))
			#y1 = centering.nona(y1)
			lims = range(y0, y1, y1, y2, y3, y4, na.rm = TRUE)
			plot(c(0:15)*3,t(y0),ylim=lims,type='b',lty=1, lwd=2.0,col='gray40', main=paste(gene,"\npvals=", pval, ", phases=", phase))
			points(c(0:15)*3,y1,type='b',col='darkblue',lty=1, lwd=2.0)
			points(c(0:3)*6,y2,type='b',col='darkred',lty=1, lwd=2.0)
			points(c(0:3)*6,y3,type='b',col='darkgreen',lty=1, lwd=2.0)
			points(c(0:3)*6,y4,type='b',col='magenta3',lty=1, lwd=2.0)
			
			abline(h=0,col='gray',lwd=3.0)
			legend('topright', legend = c('WT','Phos.WT', 'Phos.CryKO','Phos.BmalWT', 'Phos.BmalKO') , fill = c('gray40', 'darkblue','darkred','darkgreen', 'magenta3'), border = NA, bty = 'n')
		}
	}
}
dev.off()

kk = which(is.na(phos.nuclear$mapping.gene)==TRUE & phos.nuclear$gene!='')

######
###### Pipline of kinase predictions
######
data.phos = as.matrix(phos.nuclear[,c(1:16)])
data.nuclear = as.matrix(phos.nuclear[, grep('nuclear', colnames(phos.nuclear))])
ratio = data.phos - data.nuclear
res = t(apply(ratio, 1, f24_R2_alt2, t=c(0:15)*3))
res = cbind(res, qv=qvals(res[,6]))

index.sel = which(res[,1]>=6 & res[,6]<0.05 & phos.nuclear$Fasta.headers!='')
output= res[sel,]
colnames(output) = phos.nuclear$Fasta.headers[sel]

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

################# relationship between phosphorylated proteins and total proteins

gene = unique(phos.prot[,47])
#gene = gene[-which(gene=='')]
keep = c()
nb.rhy = c()
nb.rhy.total = c()
stat.mrna = c()

for(n in 1:length(gene))
{
	kk = which(phos.prot.mrna[,47]==gene[n])
	nb.rhy = c(nb.rhy, length(which(phos.prot.mrna[kk,39]<0.25)))
	nb.rhy.total = c(nb.rhy.total, length(which(phos.prot.mrna[kk,46]<0.25)))
	
	if(length(unique(phos.prot[kk,46]))!=1) {print('error');print(n);stat.mrna = rbind(stat.mrna, c(phos.prot.mrna[kk[2], c(83,86)]));}
	if(length(unique(phos.prot[kk,46]))==1) stat.mrna = rbind(stat.mrna, c(phos.prot.mrna[kk[1], c(83,86)]));
	keep = c(keep, length(kk))
}
nb.norhy = keep - nb.rhy

group = cbind(keep, nb.rhy, nb.norhy, nb.rhy.total, stat.mrna, gene)


colnames(group) = c('nb.sites', 'nb.rhythmic.sites','nb.nonrhythmic.sites', 'nb.rhymic.total.prot','mrna.rel.amp', 'mrna.Qvals', 'Gene.names')

model1 = which(group[,4]==0 & group[,2]==0)
model2 = which(group[,4]>0 & group[,2]==0)
model3 = which(group[,4]>0 & group[,2]>0)
model4 = which(group[,4]==0 & group[,2]>0)
length(model1); length(model2);length(model3); length(model4); length(c(model1,model2,model3,model4))
length(model1)/length(c(model1,model2,model3,model4));length(model2)/length(c(model1,model2,model3,model4));length(model3)/length(c(model1,model2,model3,model4));length(model4)/length(c(model1,model2,model3,model4));
cutoff.mrna = 0.05

length(which(group[model1,6]<cutoff.mrna & group[model1,5]>=0.05));
length(which(group[model2,6]<cutoff.mrna & group[model2,5]>=0.05));
length(which(group[model3,6]<cutoff.mrna & group[model3,5]>=0.05));
length(which(group[model4,6]<cutoff.mrna & group[model4,5]>=0.05));

###### add class number into the table
class = rep(0, nrow(phos.prot.mrna))
class[which(is.na(match(phos.prot.mrna[,47], gene[model1]))==F)]=1
class[which(is.na(match(phos.prot.mrna[,47], gene[model2]))==F)]=2
class[which(is.na(match(phos.prot.mrna[,47], gene[model3]))==F)]=3
class[which(is.na(match(phos.prot.mrna[,47], gene[model4]))==F)]=4
phos.prot.mrna = cbind(phos.prot.mrna,class)
#write.table(phos.prot.mrna,'Table_Phospho_Total_Protein_mRNA_Seq_all.txt',sep='\t',quote=T,col.names=T,row.name=F)

which(group[model3,6]<cutoff.mrna & group[model3,5]>=0.05)
group = as.matrix(group)
group[model3[which(group[model3,6]<cutoff.mrna & group[model3,5]>=0.05)],7]

group[model4[which(group[model4,6]<cutoff.mrna & group[model4,5]>=0.05)],7]

kk = which(phos.prot.mrna[,39]<0.25 & phos.prot.mrna[,46]<0.25 & phos.prot.mrna[,86]<0.05 & phos.prot.mrna[,83]>=0.05)
phos.prot.mrna[kk,c(33:55, 80:87)]

kk = which(phos.prot.mrna[,39]<0.25 & phos.prot.mrna[,46]>=0.25 & phos.prot.mrna[,86]<0.05 & phos.prot.mrna[,83]>=0.05)
phos.prot.mrna[kk,c(33:55, 80:87)]

############################################# group G3 genes compare the phases of total and phosphorylated proteins. 

gene.model3 = gene[model3]
newkeep = c()
for(n in 1:length(gene.model3))
{
	kk = which(phos.prot.mrna[,47]==gene.model3[n])
	print(length(kk))
	newkeep = c(newkeep, kk[which(phos.prot.mrna[kk,39]<0.25)])
	
}

phase.phos = phos.prot[newkeep, 37]
phase.prot = phos.prot[newkeep, 44]
phase.diff = phase.phos - phase.prot

for(kk in 1:length(phase.diff))
{
	test = phase.diff[kk]
	if(test>=18) phase.diff[kk] = phase.diff[kk]-24.0
	if(test<(-6.0)) phase.diff[kk] = phase.diff[kk]+24.0
}

pdf('Figure_XYX1.pdf', width=1.5, height=1.5)
#par(pty='s')
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
hist(phase.diff,breaks=c(-6:18),xlab='phase differences',ylab='number of phospho sites', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='delay [hr]',ylab='number of phospho sites', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(-6,18,by=2),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,10, by=2),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()

######## Group G4 genes comparison of phosphorylated proteins and mrna.

kk = which(phos.prot.mrna[,39]<0.25 & phos.prot.mrna[,46]>=0.25 & phos.prot.mrna[,86]<0.05 & phos.prot.mrna[,83]>=0.05)
phase.phos = phos.prot.mrna[kk,37]
phase.mrna = phos.prot.mrna[kk,84]

#gene.model2 = gene[model2]
#kk = match(gene.model2, prot.mrna[,29])
#phase.model2 = phos.prot[newkeep, 44]

############################################# rhythmic phosphorylation sites of the same genes always have the different phase

gene.all = unique(aa[,47])
ii = which(gene.all=='')
gene.all = gene.all[-ii]
phase.keep = c()
gene.keep = c()
n = which(gene.all=='Cobll1')
for(n in 1:length(gene.all))
{
	kk = which(aa[,47]==gene.all[n])
	nb.rhy.peptides = length(which(aa[kk,39]<0.25 & aa[kk,33]>=8))
	if(nb.rhy.peptides>1)
	{
		kk = kk[which(aa[kk,39]<0.25 & aa[kk,33]>=8)]
		pphase = aa[kk,37]
		phase.keep = c(phase.keep, (pphase-mean(pphase)))
		#print(pphase-mean(pphase))
		#print(aa[kk,c(33:39)])
		if(any(abs(pphase-mean(pphase)>=3))) {gene.keep = c(gene.keep, n);print(pphase);}
		#keep = c(keep, length(kk))
	}
}


pdf('Figure_SXY2.pdf', width=1.5, height=1.5)
#par(pty='s')

par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
hist(phase.keep,breaks=10,xlab='phase differeces',ylab='number of phospho sites', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(-15,15,by=5),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,50, by=10),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()




############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################









pdf('Figure_SXY1.pdf', width=1.5, height=1.5)
#par(pty='s')
ii = which(aa[,39]<0.25)

par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
hist(aa[ii,37],breaks=c(0:8)*3,xlab='phase',ylab='Number of proteins', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(0,24,by=3),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,60, by=10),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()

pdf('Figure_Lab2.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)

gene.name = 'Per2'
ii = which(aa[,47]==gene.name)
ii = ii[1]
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])

matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=3), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.,2.0, by=0.5),las=1,cex.axis = 1.0)

box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()


pdf('Figure_Lab12.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)

gene.name = 'Slc7a2'
ii = which(aa[,47]==gene.name)
ii = ii[which(aa[ii,39]<0.25)]
#ii = ii[4]
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])

matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=3), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.6,1.2, by=0.2),las=1,cex.axis = 1.0)

box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()





#############################################
#############################################
############################################# plots for some preliminary results
#############################################
#############################################
ii = which(aa[,39]<0.25)
hist(aa[ii,35],breaks=c(0:17),xlab='Peak-trough amp [log2]',ylab='Number of proteins', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(0,17,by=2),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,250, by=50),las=1,cex.axis = 1.0)

gene.name = 'Bckdha'
jj = which(aa[,47]==gene.name)
ii = jj[6] 

pdf('Figure_X3.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h',', Qval=', round(aa[ii,39], d=3), ', position::', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.7,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,10, by=2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()


#phos.prot = read.table('phospho_WT_Ref_chosen.txt',sep='\t',header=T,as.is=c(47:49,51:55))


#####################################
#####################################  Examples of phosphoproteome
#####################################

pdf('Figure_X1.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)

gene.name = 'Per2'
ii = which(aa[,47]==gene.name)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])

matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=3), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)

box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

o1 = order(aa[,39])
aa = aa[o1,]

length(which(aa[,39]<0.25))
gene.name = aa[1,47]
ii = 1 

pdf('Figure_X2.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

gene.name = aa[2,47]
ii = 2 

pdf('Figure_X3.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

gene.name = aa[3,47]
ii = 3 

pdf('Figure_X4.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

gene.name = aa[289,47]
ii = 289

pdf('Figure_X5.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3)), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.8,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()


#####################################
#####################################  Examples of proteins with different phosphoraltion dynamics 
#####################################

gene.un = unique(aa[,47])
gene.un = gene.un[-which(gene.un=='')]
keep = c()
nb.rhy = c()
for(n in 1:length(gene.un))
{
	kk = which(aa[,47]==gene.un[n])
	nb.rhy = c(nb.rhy, length(which(aa[kk,39]<0.25)))
	keep = c(keep, length(kk))
}
nb.norhy = keep - nb.rhy


n = 5 
ii = which(aa[,47]==gene.un[n])
length(kk)
gene.name = aa[ii,47]

pdf('Figure_XX2.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3), ', Position:', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.5,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

n = 204
ii = which(aa[,47]==gene.un[n])
length(kk)
gene.name = aa[ii,47]
pdf('Figure_XX3.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', pval=', round(aa[ii,38],d=4), ', Qval=', round(aa[ii,39], d=3), ', Position:', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.5,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()


n = 4
ii = which(aa[,47]==gene.un[n])
length(ii)
ii = ii[1]
gene.name = aa[ii,47]
pdf('Figure_XX4.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', Qval=', round(aa[ii,39], d=3), ', Position:', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.7,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

n = 202
ii = which(aa[,47]==gene.un[n])
length(ii)
ii = ii[3]
gene.name = aa[ii,47]
pdf('Figure_XX7.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', Qval=', round(aa[ii,39], d=3), ', Position:', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.7,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()

n = 3
ii = which(aa[,47]==gene.un[n])
length(ii)
ii = ii[3]
gene.name = aa[ii,47]
pdf('Figure_XY1.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
y=t(aa[ii,1:16])
y = 2^y
y = y/mean(y[which(!is.na(y)==T)])
lims =range(y[which(!is.na(y)==T)])
matplot(3*c(0:15), y, type='b', lty=1, pch = 15, lwd=2.0, col="darkred", ylim=lims,main=paste(gene.name,', Phase=',round(aa[ii,37],d=1), sep='','h', ', Qval=', round(aa[ii,39], d=3), ', Position:', aa[ii,50]), ylab='Phosphoprotein abundance (normalized)', cex.lab=1.0,cex.main=0.7,xlab=NA,xlim=c(0,48),axes = FALSE)
axis(1,at=3*c(0:15),cex.axis =1.0)
axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0.5,1.5, by=0.2),las=1,cex.axis = 1.0)
box()
abline(h=1,lty=2,lwd=2.0, col="darkgray")
dev.off()



myhist <- function(x, ..., breaks="Sturges",
main = paste("Histogram of", xname),
xlab = xname,
ylab = "Frequency") {
	xname = paste(deparse(substitute(x), 500), collapse="\n")
	h = hist(x, breaks=breaks, plot=FALSE)
	plot(h$breaks, c(NA,h$counts), type='S', main=main,
		 xlab=xlab, ylab=ylab, axes=FALSE, ...)
	axis(1)
	axis(2)
	lines(h$breaks, c(h$counts,NA), type='s')
	lines(h$breaks, c(NA,h$counts), type='h')
	lines(h$breaks, c(h$counts,NA), type='h')
	lines(h$breaks, rep(0,length(h$breaks)), type='S')
	invisible(h)
}
#ii = which(prot[,24]<0.25)
#ii = which(prot[,24]<0.25)
pdf('Figure_XX1.pdf', width=1.5, height=1.5)
#par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
#par(pty='s')
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
myhist((keep),breaks=c(0:40),xlab='Nb of phosphorylated sites',ylab='Number of proteins', main=NA,log='y',cex.lab=0.7,col='gray30')
dev.off()


pdf('Figure_XY2.pdf', width=1.5, height=1.5)
#par(pty='s')
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
hist(aa[,38],breaks=c(0:10)/10,xlab='Pval',ylab='Number of proteins', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(0,1,by=0.2),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,500, by=100),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()


qq = c(0:100)/100
nb.rhythmic = c()
for(n in 1:length(qq))
{
	cutoff.qq = qq[n]
	nb.rhythmic = c(nb.rhythmic,length(which(aa[,39]<cutoff.qq)))
	
}

pdf('Figure_S1F.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
plot(qq, nb.rhythmic+0.01, type='l', lty=1, lwd=2, main='', log='y',ylim=c(1,5000), xlab='FDR', ylab='Number of rhythmic peptides', col=c('blue'), cex.axis=1.0,cex.lab=1.0,cex.main=1.5)
abline(v=0.25,col='darkgray',lwd=2,lty=2)
abline(v=0.2,col='darkgreen',lwd=2,lty=2)
abline(v=0.1,col='darkgreen',lwd=2,lty=2)
dev.off()


pdf('Figure_XY4.pdf', width=2.5, height=2.5)
#par(pty='s')
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
ii = which(aa[,39]<0.25)
hist(aa[ii,37],breaks=c(0:8)*3,xlab='Phase [h]',ylab='Number of proteins', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(0,24,by=3),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,75, by=25),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()



pdf('Figure_XY3.pdf', width=2.5, height=2.5)
#par(pty='s')
par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
ii = which(aa[,39]<0.25)
hist(aa[ii,35],breaks=c(0:17),xlab='Peak-trough amp [log2]',ylab='Number of proteins', main=NA,col='gray',cex.lab=1.0,axes=FALSE)
#hist(pp[,3],breaks=c(-6:6)*2,xlab='Delay [hr]',ylab='Number of proteins', main=NA,col='gray50',cex.lab=1.0,ylim=c(0,30),xlim=c(-12,12),axes=FALSE)
axis(1,at=seq(0,17,by=2),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,250, by=50),las=1,cex.axis = 1.0)
#axis(3,lwd = 2.0)

dev.off()

gene.un = unique(aa[,47])
gene.un = gene.un[-which(gene.un=='')]
keep = c()
nb.rhy = c()
for(n in 1:length(gene.un))
{
	kk = which(aa[,47]==gene.un[n])
	nb.rhy = c(nb.rhy, length(which(aa[kk,39]<0.25)))
	keep = c(keep, length(kk))
}
nb.norhy = keep - nb.rhy

nb.rhy.phos = cbind(keep, nb.rhy, nb.norhy)
nb.rhy = nb.rhy.phos
colnames(nb.rhy) = c('nb.peptides', 'nb.rhythmic.peptides','nb.nonrhythmic.peptides')

perct.mrna = 0.15
perct.total = 0.05
perct.peptides = length(which(aa[,39]<0.25))/nrow(aa)
perct.phos1 = length(which(nb.rhy[,1]==1 & nb.rhy[,2]==1))/(length(which(nb.rhy[,1]==1)))
perct.phos2 = length(which(nb.rhy[,1]>1 & nb.rhy[,3]==0))/(length(which(nb.rhy[,1]>1)))
perct.phos3 = length(which(nb.rhy[,1]>1 & nb.rhy[,2]>0))/(length(which(nb.rhy[,1]>1)))
perct.phos.all = (length(which(nb.rhy[,2]>0)))/nrow(nb.rhy)




pdf('Figure_XY4.pdf', width=3.5, height=3.5)
#par(pty='s')
par(cex = 0.6, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)

counts = cbind(perct.mrna, perct.total, perct.peptides, perct.phos.all, perct.phos1, perct.phos2, perct.phos3)
#counts = rbind(c(1:7),counts)
barplot(t(counts), ylim=c(0,0.4),main=NA, xlab="Percentages of rhythmic features", col=c("gray30","darkblue","red1","darkred","gold4","blueviolet","darkorange"),legend = c('mRNA','total.protein','phospho.peptides','phospho.protein'), names.arg=c(1:7),beside=TRUE,axes=FALSE)
abline(0.15,0,col='gray',lwd=2.5,lty=2)
abline(0.05,0,col='gray50',lwd=2.5,lty=2)
#axis(1,at=seq(0,16,by=2),cex.axis =1.0)
#axis(1,at=49,'ZT[hr]',tick=FALSE,cex.axis =1.0)
axis(2,at = seq(0,0.3, by=0.05),las=0,cex.axis = 1.0)

dev.off()


#####################################
#####################################  relationship between phospho proteins and total proteins
#####################################
phos.prot = read.table('phospho_WT_Ref_chosen.txt',sep='\t',header=T,as.is=c(47:49,51:55))
dim(phos.prot)

#phos.prot = phos.prot[-1234,]
gene = unique(phos.prot[,47])
#gene = gene[-which(gene=='')]
keep = c()
nb.rhy = c()
nb.rhy.total = c()
class = c()
for(n in 1:length(gene))
{
	kk = which(phos.prot[,47]==gene[n])
	nb.rhy = c(nb.rhy, length(which(phos.prot[kk,39]<0.25)))
	nb.rhy.total = c(nb.rhy.total, length(which(phos.prot[kk,46]<0.25)))
	if(length(unique(phos.prot[kk,46]))!=1) {print('error');print(n);class = c(class, phos.prot[kk[2],57]);}
	if(length(unique(phos.prot[kk,46]))==1) class = c(class, phos.prot[kk[1],57])
	keep = c(keep, length(kk))
}
nb.norhy = keep - nb.rhy

nb.phos.prot = cbind(keep, nb.rhy, nb.norhy, nb.rhy.total,class)
colnames(nb.phos.prot) = c('nb.peptides', 'nb.rhythmic.peptides','nb.nonrhythmic.peptides', 'nb.rhymic.total','Class')
model1 = which(nb.phos.prot[,4]==0 & nb.phos.prot[,2]==0)
model2 = which(nb.phos.prot[,4]>0 & nb.phos.prot[,2]==0)
model3 = which(nb.phos.prot[,4]>0 & nb.phos.prot[,2]>0)
model4 = which(nb.phos.prot[,4]==0 & nb.phos.prot[,2]>0)

#################################################################### plot the results
bb = phos.prot
pdf('profiles_phospho_prot_mapping.pdf',width=12,height=8)
for(n in 1:nrow(phos.prot))
#for(n in 1:100)
{
		
	name = bb[n,48]
	y0 =bb[n,c(1:16)]
	#y0 = y0 - mean(y0[!is.na(y0)])
	x0 = bb[n,c(17:32)]
	#x0 = x0 - mean(x0[!is.na(x0)])
	#x1 = bb[n,c(33:48)]
	#x1 = x1 - mean(x1[!is.na(x1)])
	lims = range(y0[!is.na(y0)],x0[!is.na(x0)])
	plot(c(0:15)*3,y0,ylim=lims,col='darkred',type='b',lty=1,main=paste(name,", phospho.qval=", signif(bb[n,39],d=3),",  prot.qavl=", signif(bb[n,46],d=3),", position==", bb[n,50]))
	points(c(0:15)*3,x0,type='b',col='darkblue',lty=1)
	#points(c(0:15)*3,x1,type='b',col='darkblue',lty=2)
	#abline(h=0,col='gray')

}
dev.off()

