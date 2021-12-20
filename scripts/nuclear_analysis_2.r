source('f24_modified_1.0.r')
nuclear = read.table('nuclear_proteins_L_H_log2_all.txt',sep='\t',header=T)
tf = read.table('candidates_transcription_regulators_ontology.txt',sep='\t',header=T)
bb12 = read.table('proteins_with_rhythms_12h.txt',sep='\t',header=T)
mrna = read.table('RPKM_Liver_f24.txt',sep='\t',header=T,row.names=1)
#mrna = read.table('mRNA_RNA_Seq_all.txt',sep='\t', header=T, rownames=1)
#sample1 = read.table('unorm_Liver_1.txt',sep='\t', header=T,row.names=1)
#sample2 = read.table('unorm_Liver_2.txt',sep='\t',header=T,row.names=1)
#gene.ll = read.table('mean_length_gene.txt',sep='\t',header=T)

##### processing RNA-seq data (Here we choose expressed gene with a threshold of log2(RPKM)>=1 for at least 12 time points)

mrna = mrna[,c(1:24)]
colnames(mrna)[1:24] = c("1_ZT00","1_ZT02","1_ZT04","1_ZT06","1_ZT08","1_ZT10","1_ZT12","1_ZT14","1_ZT16","1_ZT18","1_ZT20","1_ZT22","2_ZT00","2_ZT02","2_ZT04","2_ZT06","2_ZT08","2_ZT10","2_ZT12","2_ZT14","2_ZT16","2_ZT18","2_ZT20","2_ZT22")
#mrna = log2(mrna)

fcn.expressed = function(x)
{
	length(which(log2(x+0.0001)>1))
}

keep = apply(mrna,1, fcn.expressed)
ii = which(keep>=12)
mrna = mrna[ii,]
mrna = as.matrix(log2(mrna+0.001))
mrna[which(rownames(mrna)=='Dbp'), ]
res.mrna = t(apply(mrna, 1, f24_R2_alt2, t=c(0:23)*2))
res.mrna = cbind(res.mrna, Q.RNA.Seq=qvals(res.mrna[,6]))
length(which(res.mrna[,7]<0.05 & res.mrna[,4]>=0.05))/nrow(mrna) ## 2 criterion for the rhythmic mrna : FDR<0.05 and rel.amplitude>=0.05 
mrna = cbind(mrna, res.mrna)
write.table(mrna,'mRNA_RNA_Seq_all_threshold_1.txt',sep='\t',quote=F,col.names=T,row.name=T)

kk = match(nuclear[,20], rownames(mrna))
length(which(is.na(kk)==T))

match.keep = c()

for(n in 1:nrow(nuclear))
{
	gene.name = nuclear[n,20]
	kk = which(rownames(mrna)==gene.name)
	if(length(kk)>1) print(n);
	if(length(kk)==1) {match.keep = c(match.keep, kk);}
	if(length(kk)<1)
	{
		#print(n)
		gene.name = unlist(strsplit(as.character(gene.name), ";"))
		keep = c()
		for(name.ii in gene.name)
		{
			keep = c(keep, which(rownames(mrna)==name.ii))
		}
		keep = unique(keep)
		if(length(keep)==0) match.keep = c(match.keep, NA)
		if(length(keep)>1) 
		{
			#print(n)
			match.keep = c(match.keep, keep[which(mrna[keep,31]==min(mrna[keep,31]))])
			if(length(keep[which(mrna[keep,31]==min(mrna[keep,31]))])>1) print(length(keep))
		}
		if(length(keep)==1) match.keep = c(match.keep, keep)
	}
}

xx = cbind(nuclear, mrna[match.keep,], rownames(mrna[match.keep,]))
#Q.RNA.Seq = res.mrna[match.keep,7]
#phos.prot = cbind(phos.prot, Q.RNA.Seq)
colnames(xx)[60]='RNA.Gene.names'
nuclear.mrna = xx

class = rep(0, nrow(nuclear))
model1 = which(nuclear.mrna[,28]>=0.05 & (nuclear.mrna[,59]>=0.05|nuclear.mrna[,56]<=0.05))
model2 = which(nuclear.mrna[,28]>=0.05 & (nuclear.mrna[,59]<0.05 & nuclear.mrna[,56]>0.05))
model3 = which(nuclear.mrna[,28]<0.05 & (nuclear.mrna[,59]<0.05 & nuclear.mrna[,56]>0.05))
model4 = which(nuclear.mrna[,28]<0.05 & (nuclear.mrna[,59]>=0.05|nuclear.mrna[,56]<=0.05))

class[model1] = 1
class[model2] = 2
class[model3] = 3
class[model4] = 4


nuclear.mrna = cbind(nuclear.mrna, classes=class)
write.table(nuclear.mrna,'Table_nuclear_Protein_mRNA_Seq_all.txt',sep='\t',quote=F,col.names=T,row.name=F)






output.pol2 = read.table('Pol2_TSS.txt', sep='', header=T,row.names=1)
annot.ens = read.table('transcript_mapping.txt', sep='\t')
output.pol2 = as.matrix(output.pol2)                                ### Pol2 output and TF motifs
output.pol2 = log2(output.pol2+1)
res.pol2 = t(apply(output.pol2, 1, f24_R2_alt2, t=(c(0:6)*4+2)))
res.pol2 = cbind(res.pol2, Qval=qvals(res.pol2[,6]))
#### here we need to filter the background
cutoff = 7.0
ii = which(res.pol2[,2]>cutoff)
pol2 = output.pol2[ii,]
res.pol2 = res.pol2[ii,]

mm = match(rownames(res.pol2), annot.ens[,1])
length(which(is.na(mm)==T))
res.pol2 = res.pol2[which(!is.na(mm)==T),]
pol2 = pol2[which(!is.na(mm)==T),]
rownames(res.pol2) = annot.ens[mm[which(!is.na(mm)==T)],3]
rownames(pol2) = annot.ens[mm[which(!is.na(mm)==T)],3]
pol2 = cbind(pol2, res.pol2)


prot0 = read.table('table_all_proteins_nb_timepoints_L_H_log2.txt',sep="\t",header=T,as.is = c(26:29))  ###### total proteins measured with MS which is considered as cytoplasmic proteins
prot = prot0
jj = match(nuclear$Gene.names, prot[,29])
length(which(!is.na(jj)==T))

prot.mapping = prot[jj,]
nuclear.mapping = cbind(nuclear.mrna, prot.mapping)

jj = match(nuclear$Gene.names, rownames(pol2))
length(which(!is.na(jj)==T))
pol2.mapping = pol2[jj,]

nuclear.mapping = cbind(nuclear.mapping, pol2.mapping)

ii = grep('Eif', rownames(nuclear.mapping))

######################### check the relationship between nucelar proteins and pol2, mrna and total proteins, probably also nascent mrna

ii = which(nuclear.mapping[,103]<0.05 & nuclear.mapping[,28]<0.05)

circular.error=function(a,x,y)
{
	sum(1-cos(2*pi/24*(y-x-a)))
}

fit=nlm(circular.error, 4, x=nuclear.mapping[ii,82], y=nuclear.mapping[ii,26])
plot(nuclear.mapping[ii,c(82,26)])
abline(fit$estimate-48,1,col='red',lwd=2.0)
abline(fit$estimate-24,1,col='red',lwd=2.0)
abline(0,1,col='blue',lwd=2.0)

ii = which(nuclear.mapping[,103]<0.05 & nuclear.mapping[,28]<0.05) ### phase of pol2 vs. nuclear proteins
jj = which(nuclear.mapping[,59]<0.05 & nuclear.mapping[,28]<0.05) ### phase of mrna vs nuclear proteins




#kk = match(rownames(sample1), rownames(sample2))
#a = cbind(sample1[kk,], sample2)
#colnames(a)[1:24] = c("1_ZT00","1_ZT02","1_ZT04","1_ZT06","1_ZT08","1_ZT10","1_ZT12","1_ZT14","1_ZT16","1_ZT18","1_ZT20","1_ZT22","2_ZT00","2_ZT02","2_ZT04","2_ZT06","2_ZT08","2_ZT10","2_ZT12","2_ZT14","2_ZT16","2_ZT18","2_ZT20","2_ZT22")
#table = a[,c(1:24)]
#table = as.matrix(table)
#find_expressed = function(x) 
#{
#	length(which(x>=5))
#}
#nb.test = apply(table, 1, find_expressed)
#keep = which(nb.test>=24) #### here we select the genes in which reads in all time points >=5 
#table = table[keep,]
#jj = match(rownames(table),gene.ll[,1])
#gene.ll = gene.ll[jj,]
#genes.name = a[keep, c(13,14)]

#ii = match(rownames(table), rownames(mrna))
mrna = mrna[ii,]
mrna = mrna[,c(1:24)]
colnames(mrna)[1:24] = c("1_ZT00","1_ZT02","1_ZT04","1_ZT06","1_ZT08","1_ZT10","1_ZT12","1_ZT14","1_ZT16","1_ZT18","1_ZT20","1_ZT22","2_ZT00","2_ZT02","2_ZT04","2_ZT06","2_ZT08","2_ZT10","2_ZT12","2_ZT14","2_ZT16","2_ZT18","2_ZT20","2_ZT22")
mrna = log2(mrna)
res = t(apply(as.matrix(mrna[,c(1:24)]), 1, f24_R2_alt2, t=c(0:23)*2))
qv = qvals(res[,6])
res = cbind(res,qv)
mrna = cbind(mrna,res)

mapping = c()
#trace = c()
keep = c()
for(n in 1:nrow(nuclear))
{
	gene.name =  unlist(strsplit(as.character(nuclear[n,20]), ";"))
	gene.name = unique(gene.name)
	if(length(gene.name)<1) { print("NO name");mapping = c(mapping, NA);}
	if(length(gene.name)>=1)
	{
		trace = c()
		for(kk in 1:length(gene.name))
		{
			trace = c(trace, which(rownames(mrna)==gene.name[kk]))
		}
		trace = unique(trace)
		if(length(trace)==0) {mapping = c(mapping, NA); keep = c(keep, n)}
		if(length(trace)==1) mapping = c(mapping, trace)
		if(length(trace)>1) mapping = c(mapping, trace[which(mrna[trace,31]==min(mrna[trace,31]))])
	}
	
}

nuclear.mrna = cbind(nuclear, mrna[mapping,], rownames(mrna)[mapping])
colnames(nuclear.mrna)[c(22:28, 53:60)] = c('nb.timepoints.nuclear','mean.nuclear','amp.nuclear','relamp.nuclear','phase.nuclear','pval.nuclear','qval.nuclear',
'nb.timepoints.mrna','mean.mrna','amp.mrna','relamp.mrna','phase.mrna','pval.mrna','qval.mrna','Gene.names.mRNA')
plot(nuclear.mrna[,c(58,27)],log='xy',xlim=c(0.00001,1),ylim=c(0.00001,1),ylab='pval of nuclear proteins', xlab='pval of mrna') 
abline(h=0.01,col='green',lwd=2.5)
abline(v=0.01,col='green',lwd=2.5)
abline(0,1,col='red',lwd=2.5)

model1 = which(nuclear.mrna[,59]>=0.05 & nuclear.mrna[,28]>=0.05)
model2 = which(nuclear.mrna[,59]<0.05 & nuclear.mrna[,28]>=0.05)
model3 = which(nuclear.mrna[,59]<0.05 & nuclear.mrna[,28]<0.05)
model4 = which(nuclear.mrna[,59]>=0.05 & nuclear.mrna[,28]<0.05)

all = c(model1, model2, model3, model4)

length(model1)/length(all);length(model2)/length(all);length(model3)/length(all);length(model4)/length(all);

#write.table(nuclear.mrna,'nuclear_proteins_mrna.txt',sep='\t',quote=F,col.names=T,row.names=F)


###############################################
############################################### gene ontology analysis
###############################################

id.all = c()

for(n in 1:nrow(nuclear))
{
	test = nuclear[n,18]
	test = unlist(strsplit(as.character(test), ";"))
	id.all = c(id.all, test)
	#lines = c(rownames(trans)[n], exons.x[k,8], exons.x[k,8]+(overlap[ii[nnb],3]-start), gene.name,0, trans[n,7])
	#redregion = rbind(redregion, lines)
	#lines = paste(lines,collapse='\t')
	#write(lines,file='output.txt',append=TRUE)
	
}
write.table(id.all,'protein_id_all_uniprot.txt',sep='\t',quote=F,col.names=F,row.names=F)

id.rhythmic = c()
for(n in c(model3, model4))
{
	test = nuclear[n,18]
	test = unlist(strsplit(as.character(test), ";"))
	id.rhythmic = c(id.rhythmic, test)
}
write.table(id.rhythmic,'protein_id_rhythmic_uniprot.txt',sep='\t',quote=F,col.names=F,row.names=F)

id.model3 = c()
for(n in c(model3))
{
	test = nuclear[n,18]
	test = unlist(strsplit(as.character(test), ";"))
	id.model3 = c(id.model3, test)
}
write.table(id.model3,'protein_id_model3_uniprot.txt',sep='\t',quote=F,col.names=F,row.names=F)

id.model4 = c()
for(n in c(model4))
{
	test = nuclear[n,18]
	test = unlist(strsplit(as.character(test), ";"))
	id.model4 = c(id.model4, test)
}
write.table(id.model4,'protein_id_model4_uniprot.txt',sep='\t',quote=F,col.names=F,row.names=F)



###############################################
############################################### rhythmic proteins with different functions
###############################################
tf = read.table('candidates_transcription_regulators_ontology.txt',sep='\t',header=T)
kinase = read.csv('kinase.csv',header=F)
phosphatase = read.csv('phosphatases.csv',header=F);
rna.processing = read.csv('genes_RNA_processing.csv',header=F)
chromatin = read.csv('genes_chromatin_regulators.csv', header=F)
histone = read.csv('genes_histones_motification.csv',header=F)
splicesome = read.csv('genes_splicesomes.csv',header=F)
translation.factors = read.csv('translaton_factors.csv',header=F)

tf = unique(tf[,3])
kk = match(tf, nuclear.mrna[,60])
kk = kk[which(!is.na(kk)==T)]
length(tf);length(kk);length(which(nuclear.mrna[kk,28]<0.05));length(which(nuclear.mrna[kk,28]<0.05))/length(kk)

regulator = translation.factors
regulator = unique(regulator[,3])
kk = match(regulator, nuclear.mrna[,60])
kk = kk[which(!is.na(kk)==T)]
length(regulator);length(kk);length(which(nuclear.mrna[kk,28]<0.05));length(which(nuclear.mrna[kk,28]<0.05))/length(kk)

ii = grep('Eif', rownames(nuclear.mapping))
kk = unique(c(ii,kk))
kk = kk[order(kk)]

kk = c(kk, grep('40S', nuclear.mapping[,19]), grep('60S', nuclear.mapping[,19]))
kk = unique(kk)
cutoff = 0.05

which(nuclear.mapping[kk,28]<cutoff)
length(which(nuclear.mapping[kk,28]<cutoff))
ii = kk[which(nuclear.mapping[kk,28]<cutoff)]
nuclear.mapping[ii,c(22:28)]
hist(nuclear.mapping[ii,26],breaks=c(0:8)*3)

jj = c(1:530)




###### plots of proteins with different functions 
regulator = tf
regulator = unique(regulator[,3])
kk = match(regulator, nuclear.mrna[,60])
kk = kk[which(!is.na(kk)==T)]
subgroup = nuclear.mrna[kk,]
o2 = order(subgroup$qval.nuclear)
subgroup = subgroup[o2,]

pdf('Nulcear_tf.pdf',width=10,height=8)
for(n in 1:nrow(subgroup))
#for(n in 1:20)
{
	if(subgroup[n,22]>3)
	{
		name = subgroup[n,60]
		y0 = as.numeric(subgroup[n,c(1:16)])
		y0 = y0-mean(y0[which(!is.na(y0)==T)])
		time = c(0:15)*3
		yy0 = as.numeric(subgroup[n, c(29:52)])
		yy0 = yy0-mean(yy0)
		lims = range(c(y0[which(!is.na(y0)==T)], yy0))
		plot(time,y0,ylim=lims,col='darkred',type='b',lwd=2.0,ylab='nuclear.prot abundance', xlab='ZT[h]',lty=1, 
			 main=paste(name,", phase.prot=", signif(subgroup[n,26],d=3),", qval.prot=", signif(subgroup[n,28],d=3), ",phase.mrna=", signif(subgroup[n,57],d=3),", qval.mrna=", signif(subgroup[n,59],d=3)))
		#plot(time,y0,ylim=lims,col='darkred',type='b',lty=1)
		lines(c(0:23)*2, yy0, col='darkblue', type='b',lwd=2.0)
		#abline(h=res.factors[n,2],col='gray',lwd=2.0)
	}
}
dev.off()




############################################### ###############################################
############################################### ###############################################
############################################### ###############################################
############################################### ###############################################
############################################### ###############################################
############################################### ###############################################
############################################### ###############################################


require(preprocessCore)
a = read.table('proteinGroups_Filtered_all_log2.txt',header=T,sep='\t',na.strings="NaN")

aa = a[,c(1:16, 37:41)]
bb = -aa[,c(1:16)]
aa = cbind(bb,aa[,c(17:21)])
colnames(aa)[1:16] = c("ZT00.WT","ZT03.WT","ZT06.WT","ZT09.WT","ZT12.WT","ZT15.WT","ZT18.WT","ZT21.WT","ZT24.WT","ZT27.WT","ZT30.WT","ZT33.WT","ZT36.WT","ZT39.WT","ZT42.WT","ZT45.WT")
source('f24_modified_1.0.r')

res = t(apply(as.matrix(aa[,c(1:16)]), 1, f24_R2_alt2, t=c(0:15)*3))
qv = qvals(res[,6])
res = cbind(res,qv)

res2 = t(apply(as.matrix(aa[,c(1:16)]), 1, f24_R2_alt2, t=c(0:15)*3, period=12))
qv2 = qvals(res2[,6])
res2 = cbind(res2,qv2)

#res3 = t(apply(as.matrix(aa[,c(1:16)]), 1, f24_R2_alt2, t=c(0:15)*3, period=6))
#qv2 = qvals(res2[,6])
#res2 = cbind(res2,qv2)

ii = order(res[,7])
res = res[ii,]
aa = aa[ii,]

aa = cbind(aa, res)
write.table(aa,'nuclear_proteins_L_H_log2_all.txt',sep='\t',quote=F,col.names=T,row.names=F)


###############################################
############################################### condidats of transcriptional factors
###############################################
tf = read.csv('Transcrition_factor_mouse.csv',header=F)
tf = tf[,-4]

ttf = tf
for(n in 1:nrow(tf))
{
	kk = which(ttf[n,]=='')
	if(length(kk)>0) 	ttf[n,kk] = NA
}

prot.tf = unique(ttf[,3])

tfb = read.csv('binding_tf.csv',header=F)

tfbb = tfb
for(n in 1:nrow(tfbb))
{
	kk = which(tfbb[n,]=='')
	if(length(kk)>0) 	tfbb[n,kk] = NA
}
tfbb =tfbb[,-4]
tfbb = tfbb[,-8]
tfbb = tfbb[,c(1,2,3,4,7:9)]

write.table(tfbb,'candidates_transcription_regulators_ontology.txt',sep='\t',quote=F,col.names=F,row.names=F)

prot.tf = unique(tfbb[,3])
prot.tf1 =  unique(ttf[,3])


kk = match(prot.tf, aa[,20])
length(kk)
length(kk[which(is.na(kk)==F)])

factors = aa[kk[which(is.na(kk)==F)],]
res.factors = res[kk[which(is.na(kk)==F)],] 

o2 = order(res.factors[,7])
factors = factors[o2,]
res.factors = res.factors[o2,]

write.table(factors,'Nuclear_transcription_regulators_ontology_detected.txt',sep='\t',quote=F,col.names=T,row.names=F)


pdf('Nulcear_transcriptional_factors.pdf',width=10,height=8)
for(n in 1:nrow(factors))
#for(n in 1:20)
{
	if(res.factors[n,1]>3)
	{
		name = factors[n,20]
		y0 = factors[n,c(1:16)]
		time = c(0:15)*3
		lims = range(y0[which(!is.na(y0)==T)])
		plot(time,y0,ylim=lims,col='darkred',type='b',lwd=2.0,ylab='protein abundance', xlab='ZT[h]',lty=1, main=paste(name," phase=", signif(res.factors[n,5],d=3),", pval=", signif(res.factors[n,6],d=3),", qval=", signif(res.factors[n,7],d=3)))
		#plot(time,y0,ylim=lims,col='darkred',type='b',lty=1)
		abline(h=res.factors[n,2],col='gray',lwd=2.0)
	}
}
dev.off()



pdf('Nulcear_proteins_profiles_all.pdf',width=12,height=8)
for(n in 1:nrow(aa))
#for(n in 1:20)
{
	if(res[n,1]>3)
	{
		name = aa[n,20]
		y0 = aa[n,c(1:16)]
		time = c(0:15)*3
		lims = range(y0[which(!is.na(y0)==T)])
		plot(time,y0,ylim=lims,col='darkred',type='b',lwd=2.0,ylab='protein abundance', xlab='ZT[h]',lty=1, main=paste(name," phase=", signif(res[n,5],d=3),", pval=", signif(res[n,6],d=3), ", qval=", signif(res[n,7],d=3)))
		#plot(time,y0,ylim=lims,col='darkred',type='b',lty=1)
		abline(h=res[n,2],col='gray',lwd=2.0)
	}
}
dev.off()

kk = which(res2[,6]<0.05)
bb = aa[kk,]
res2 = res2[kk,]
o1 = order(res2[,6])
res2 = res2[o1,]
bb = bb[o1,]
bb = cbind(bb, res2)

colnames(bb)[29:34] = c('nb.timepoints.f12','mean.f12','amp.f12','relamp.f12','phase.f12', 'pval.f12')
bb = bb[,-c(29,30,35)]

write.table(bb,'proteins_with_rhythms_12h.txt',sep='\t',quote=F,col.names=T,row.names=F)


pdf('Nulcear_proteins_profiles_12h_rhythm.pdf',width=12,height=8)
for(n in 1:nrow(bb))
#for(n in 1:20)
{
	if(res2[n,1]>3)
	{
		name = bb[n,20]
		y0 = bb[n,c(1:16)]
		time = c(0:15)*3
		lims = range(y0[which(!is.na(y0)==T)])
		plot(time,y0,ylim=lims,col='darkred',type='b',lwd=2.0,ylab='protein abundance', xlab='ZT[h]',lty=1, main=paste(name," phase=", signif(res2[n,5],d=3),", pval=", signif(res2[n,6],d=3)))
		#plot(time,y0,ylim=lims,col='darkred',type='b',lty=1)
		abline(h=res2[n,2],col='gray',lwd=2.0)
	}
}
dev.off()



qq = c(0:100)/100
nb.rhythmic = c()
for(n in 1:length(qq))
{
	cutoff.qq = qq[n]
	nb.rhythmic = c(nb.rhythmic,length(which(res[,7]<cutoff.qq)))
	
}

pdf('Figure_S1F.pdf', width=2.5, height=2.5)
par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
plot(qq, nb.rhythmic+0.01, type='l', lty=1, lwd=2, main='', log='y',ylim=c(1,5000), xlab='FDR', ylab='Number of rhythmic proteins', col=c('blue'), cex.axis=1.0,cex.lab=1.0,cex.main=1.5)
abline(v=0.05,col='darkred',lwd=2,lty=2)
abline(v=0.2,col='darkgreen',lwd=2,lty=2)
abline(v=0.1,col='darkgreen',lwd=2,lty=2)
dev.off()



ii = which(res[,7]<0.25)


list = c('Arntl','Dbp', 'Arntl2','Per1', 'Per2', 'Per3','Cry1', 'Cry2','Bhlhe41','Tef', 'Bhlhe40','Hlf','Alb', 'Nr1d1','Nr1d2','Npas2','Rhob','Rorg', 'Nfil3','Ck1d/e')
time = c(0,6,12,18,24,30,36,42,45)
pdf('first_glimpse_Nulcear_proteins_circadian_from_all.pdf',width=12,height=8)
for(n in 1:length(list))
{
	ii = which(aa[,13]==list[n])
	if(length(ii)==0) {print(list[n]); print('No found');}
	if(length(ii)>0)
	{
		y0 = aa[ii,c(1:9)]
		lims = range(y0[which(!is.na(y0)==T)])
		plot(time,y0,ylim=lims,col='darkred',type='b',lty=1, main=paste(list[n],", pval=", signif(res[ii,6],d=2), ", qval=", signif(res[ii,7],d=2), ", phase=", signif(res[ii,5],d=2), "h"))
		#plot(time,y0,ylim=lims,col='darkred',type='b',lty=1)
		abline(h=res[n,2],col='gray')
		
	}
}
dev.off()

##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################

Per1
Dbp
Rora
Tef
Cry2
ENSMUSG00000055866	1	93312561	93355907	-1	Per2
ENSMUSG00000030256	6	145806763	145814078	-1	Bhlhe41
ENSMUSG00000026077	1	39250576	39420079	1	Npas2
ENSMUSG00000040187	6	146744577	146782051	1	Arntl2
ENSMUSG00000020889	11	98629246	98636647	-1	Nr1d1
ENSMUSG00000003949	11	90197849	90252231	-1	Hlf
ENSMUSG00000029238	5	76639992	76733817	-1	Clock
ENSMUSG00000036192	19	19005095	19185686	-1	Rorb
ENSMUSG00000055116	7	120350979	120457636	1	Arntl
ENSMUSG00000056749	13	53062578	53076442	-1	Nfil3
ENSMUSG00000030103	6	108610623	108616919	1	Bhlhe40
ENSMUSG00000028957	4	150377761	150418774	-1	Per3
ENSMUSG00000020038	10	84594447	84647799	-1	Cry1
ENSMUSG00000021775	14	19036568	19071641	-1	Nr1d2
ENSMUSG00000029878	12	75398461	75401455	1	Dbpht2







aa = a[,c(42:57,38:41,34:37,30:33,5,1:4,9)]
bb = -aa[,c(1:28)]
aa = cbind(bb,aa[,c(29:34)])
colnames(aa)[1:28] = c("ZT00.WT","ZT03.WT","ZT06.WT","ZT09.WT","ZT12.WT","ZT15.WT","ZT18.WT","ZT21.WT","ZT24.WT","ZT27.WT","ZT30.WT","ZT33.WT","ZT36.WT","ZT39.WT","ZT42.WT","ZT45.WT","ZT00.Cry.KO","ZT06.Cry.KO","ZT12.Cry.KO","ZT18.Cry.KO","ZT00.Bmal.WT","ZT06.Bmal.WT","ZT12.Bmal.WT","ZT18.Bmal.WT","ZT00.Bmal.KO","ZT06.Bmal.KO","ZT12.Bmal.KO","ZT18.Bmal.KO")
#write.table(aa,'Bmal_Cry_WT_KO_all_samples.txt',sep='\t',quote=F,col.names=T)

write.table(aa,'Phospho_all.txt',sep='\t',quote=F,col.names=T)

data = as.matrix(aa[,c(1:16)])
res = c()
for(n in 1:nrow(data))
{
	test = data[n,]
	t = which(!is.na(test)==TRUE)
	test = test[t]
	xx = f24_R2_alt(test,(t-1)*3)
	res = rbind(res,xx)
}

o = order(res[,6])
res = res[o,]
aa = aa[o,]

pdf('first_glimpse.pdf',width=12,height=8)
for(n in 1:nrow(aa))
#for(n in 1:200)
{
	if(length(which(is.na(aa[n,c(1:28)])==TRUE))<7)
	{
	name = aa[n,29]
	y0 = aa[n,c(1:16)]
	y00 = aa[n,c(9:16)]
	y1 = aa[n,c(17:20)]
	y2 = aa[n,c(21:24)]
	y3 = aa[n,c(25:28)]
	y0 = y0 - mean(y0[!is.na(y0)])
	y00 = y00 - mean(y00[!is.na(y00)])
	y1 = y1 - mean(y1[!is.na(y1)])
	y2 = y2 - mean(y2[!is.na(y2)])
	y3 = y3 - mean(y3[!is.na(y3)])
	lims = range(y0[!is.na(y0)],y00[!is.na(y00)],y1[!is.na(y1)],y2[!is.na(y2)],y3[!is.na(y3)])
	plot(c(0:15)*3,y0,ylim=lims,col='darkblue',type='b',lty=1,main=paste(name,", pval=", signif(res[n,6],d=2)))
	#points(c(0:7)*3,y00,col='darkblue',type='b',lty=1)

	points(c(0:3)*6,y3,type='b',col='darkblue',lty=2)
	
	points(c(0:3)*6,y1,type='b',col='red',lty=1)
	points(c(0:3)*6,y2,type='b',col='red',lty=2)
	abline(h=0,col='gray')
	}
}
dev.off()

prot0 = read.table('table_all_proteins_nb_timepoints_L_H_log2.txt',sep="\t",header=T,as.is = c(26:29))
prot = prot0
prot.ko = read.table('Bmal_Cry_WT_KO_all_samples.txt',sep="\t",header=T,as.is=c(13:17))

nb = c()
for(n in 1:nrow(prot.ko))
{
	test = prot.ko[n,c(1:12)]
	nb = c(nb, length(which(!is.na(test)==TRUE)))
}
prot.ko = prot.ko[which(nb==12),]
jj = match(prot.ko[,16], prot[,29])
length(which(is.na(jj)==TRUE))
map = c()
lost = c()
for(n in 1:nrow(prot))
{
	ii = which(prot.ko[,16]==prot[n,29])
	prot.ko[ii,]
	map = c(map, length(ii))
	#ii = which(prot.ko[,16]==prot[n,29] & prot.ko[,14]==prot[n,27])
	#if(length(ii)>0) map = c(map, ii)	
	#if(length(ii)<1) {map = c(map,NA);lost = c(lost, n)}
	
}




a=read.table('proteinGroups.txt', sep="\t", header=T, na.strings="Non Numérique")
b= read.table('Gachon_Cry_Bmal_12samples_4times_1quant_REV_CON_filter.txt',sep="\t",header=T,na.strings="Non Num\303\251rique")

#ii=(0:15)*4+113 #"Ratio.H.L.ZT0"
aa=b[,c(1:12)]		# normalized

aa=as.matrix(-(aa))
aa = cbind(aa,b[,c(32:36)])
aa = aa[,c(5:8,1:4,9:17)]
colnames(aa)[1:12] = c("ZT00.Bmal.WT","ZT06.Bmal.WT","ZT12.Bmal.WT","ZT18.Bmal.WT","ZT00.Bmal.KO","ZT06.Bmal.KO","ZT12.Bmal.KO","ZT18.Bmal.KO","ZT00.Cry.KO","ZT06.Cry.KO","ZT12.Cry.KO","ZT18.Cry.KO")
write.table(aa,'Bmal_Cry_WT_KO_all_samples.txt',sep='\t',quote=F,col.names=T)

ii = match(aa[,16], mutant[,1])

prot0 = read.table('table_all_proteins_nb_timepoints_L_H_log2.txt',sep="\t",header=T,as.is = c(26:29))
prot = prot0
mutant = read.table('table_cry_ko_complet.txt',sep="\t",header=T,as.is=c(1:3))

ii = match(aa[,16],prot[,29])
ref = prot[ii,c(1:16)]
names = prot[ii,29]

