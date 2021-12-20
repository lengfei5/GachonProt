nuclear = read.table('Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',as.is=c(2,3), header=TRUE, sep='\t')
#load(file='Tables_DATA/Table_Nuclear_mRNA_Total_Nascent.Rdata')
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v3.txt', sep='\t', header=TRUE)
o2 = order(table.sx$qv)
table.sx = table.sx[o2,]
load(file='Tables_DATA/Table_Nuclear_total_mRNA_pol2_v2.Rdata')
nuclear.all$qv.mrna = qvals(nuclear.all$pval.mrna)
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/Annotation_Used_All_V2.Rdata')
#load(file='Tables_DATA/Table_mRNA_CRY_WT_NRF.Rdata')
source('functions_nuclear.R')
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/mRNA_WT_RF_total_RNA_seq_Cedric.Rdata')
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/List_genes_expresed_mouse_liver.Rdata')
### define genes epxressed in mouse liver
Define.genes.liver = FALSE
if(Define.genes.liver)
{
    load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/mRNA_WT_RF_total_RNA_seq_Cedric.Rdata')
    kk = which(mrna$mean.mrna>-1)
    #total.names = read.table('Tables_DATA/Table_total_proteins_index_names.txt', sep='\t', header=TRUE, as.is=c(2))
    gene.expressed = unique(c(as.character(mrna[kk, 1])))
    save(gene.expressed, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/List_genes_expresed_mouse_liver.Rdata')

}

Table.Sup = FALSE
if(Table.Sup)
{
    xx = nuclear[, c(20, 19, 17, 18, 1:16, 34:37, 22:28)]
    colnames(xx)[25:31] = paste(colnames(xx)[25:31], '.WT', sep='')
    ### compute the mean and SEM
    kk = grep('ZT', colnames(xx))
    kk = kk[c(1:16)]
    dat = as.matrix(xx[,kk])
    newdat = matrix(NA, ncol = 16, nrow = nrow(dat))
    source('functions_nuclear.R')
    for(n in c(1:8))
    {
      test = dat[, c(n, (n+8))]
      newdat[, (n*2-1)] = apply(test, 1, mean.nona);
      newdat[, (n*2)] = apply(test, 1, sme.nona);
    }
    colnames(dat)[9:16] = paste('ZT', c(0:7)*3, '.WT.Rep2', sep='')
    colnames(dat)[1:8] = paste('ZT', c(0:7)*3, '.WT.Rep1', sep='')
    newdat = data.frame(newdat)
    colnames(newdat)[(c(1:8)*2-1)] = paste('ZT', c(0:7)*3, '.WT.Mean', sep='')
    colnames(newdat)[(c(1:8)*2)] = paste('ZT', c(0:7)*3, '.WT.SEM', sep='')
    
    xx = data.frame(xx[, c(1:4)], newdat, dat, xx[21:31], stringsAsFactors = FALSE)
    write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/paper/Revision/Revision_2/Table_S1_nuclear_proteins_all_Mean_SEM.txt', 
                sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    Revision.table.S1.MaxQuant.output = FALSE
    if(Revision.table.S1.MaxQuant.output)
    {
      yy = read.csv(file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Revision/Olds/proteinGroups_NBPEPTIDE_UNIQUE_SILAC_VAR.csv', header=FALSE, 
                    stringsAsFactors=FALSE)
      colnames(yy) = (yy[1,])
      yy = yy[-1, ]
      mm = match(xx$Protein.IDs, yy$`Protein IDs`)
      length(which(is.na(mm)))
      
      #yy = data.frame(xx, yy[mm, -c(1, 2, 6, 7)])
      yy = yy[mm, ]
      yy = yy[, c(7, 6, 1:5, 8:134)]
      write.table(yy, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Revision/Table_S1_Maxquant_output_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    }
    
    Revision.correlaiton.replicates = FALSE
    if(Revision.correlaiton.replicates)
    {
      xx = nuclear[, c(20, 19, 17, 18, 1:16, 34:37, 22:28)]
      colnames(xx)[25:31] = paste(colnames(xx)[25:31], '.WT', sep='')
      
      #### correlations between replicates
      kk = which(xx$nb.timepoints.WT>=8 & xx$qv.WT<1)
      data = as.matrix(xx[kk, c(5:20)])
      rownames(data) = xx$Gene.names[kk]
      
      ### find 300 outliers and analyze their phases and amplitudes with and without ZT21
      rr0 = c()
      for(n in 1:7){rr0 = c(rr0, (data[,n]-data[, (n+8)]));}
      hist(rr0, breaks=200)
      #data = data[which(index.outliers==0), ]
      rr = data[, 8]-data[,16]
      cutoff = quantile(rr0, probs = 0.99, na.rm = TRUE)
      index.outliers = kk[which(rr>cutoff)]
      length(index.outliers)
      
      #data = data[which(rr<cutoff), ]
      data1 = as.matrix(xx[index.outliers, c(5:20)]);
      data1[,8] = NA;
      res1 = t(apply(data1, 1, f24_R2_alt2, t=c(0:15)*3))
      res1[, 4] = t(apply(2^data1, 1, f24_R2_alt2, t=c(0:15)*3))[,4]
      #res = cbind(res, qv=qvals(res[,6]))
      
      test = data.frame(xx$Gene.names[index.outliers],  xx$qv.WT[index.outliers], xx$pval.WT[index.outliers], res1[, 6], xx$phase.WT[index.outliers], res1[,5], 
                        xx$relamp.WT[index.outliers], res1[,4], xx$ZT21.WT[index.outliers], xx$ZT45.WT[index.outliers], stringsAsFactors = FALSE)
      colnames(test) = c('gene', 'qval', 'pval.init','pval.refit', 'phase.init', 'phase.refit', 'relamp.init', 'relamp.refit', 'ZT21', 'ZT45')
      
      pdf.name = paste("myplots/FIGURES/Correlation_replicates_all_timepoints.pdf", sep='')
      pdf(pdf.name, width=4.0, height=2.2)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      
      par(mfrow=c(2, 4))
      cex = 0.02;pch=1;
      lims = c(-2., 2.)
      xlab=NA; ylab=NA;
      par(mgp = c(0,0.25,0.), mar = c(1.2,0.4,1.,0.), pty='s');
      for(n in 1:7)
      {
       plot(data[, c(n, (n+8))], cex=cex, xlab=xlab, ylab=ylab, pch=pch,
             xlim=lims, ylim=lims, main=paste('ZT', (n-1)*3, sep=''), cex.main=0.7, axes = FALSE);
        jj = which(!is.na(data[,n]==TRUE & !is.na(data[, (n+8)])))
        #test = (cor.test(data[jj,n], data[jj, (n+8)], alternative='greater'))
        if(n==1|n==5) axis(2, at=seq(-2, 2, by=2), cex=0.7);
        if(n>=5) axis(1, at=seq(-2, 2, by=2), cex=0.7);
        box();
        text(-1., 1.5, cex=0.7, paste('R = ', signif(cor(data[,n], data[, (n+8)], use="na.or.complete"), d=2), sep=''), col='red')
        abline(0, 1, lwd=1., col='red')
      }
      #par(mgp = c(0,0.3,0.), mar = c(1.2,0.2,1.0,0.1), pty='s')
      plot(data[which(rr<=cutoff), c(8, 16)], cex=cex, xlab=xlab, ylab=ylab, xlim=lims, ylim=lims, main='ZT21', pch=pch,cex.main=0.7, axes=FALSE);
      points(data[which(rr>cutoff), c(8, 16)], cex=cex, col='darkorange', pch=pch);
      axis(1, at=seq(-2, 2, by=2), cex=0.7);
      box()
      text(-1, 1.5, cex=0.7, paste('R = ', signif(cor(data[ ,8], data[, 16], use="na.or.complete"), d=2), sep=''), col='red')
      abline(0, 1, lwd=1.5, col='red')
      
      dev.off()
      
      plot(data[, c(8, 5)], cex=cex, xlab='ZT21-Rep1', ylab='ZT12-Rep1', xlim=c(-4, 4), ylim=c(-4, 4), main=NA);
      text(-2, 3, paste('R = ', signif(cor(data[,8], data[, 5], use="na.or.complete"), d=2), sep=''), col='red')
      #abline(0, 1, lwd=1.5, col='red')
      plot(data[, c(16, 5)], cex=cex, xlab='ZT21-Rep2', ylab='ZT12-Rep1', xlim=c(-4, 4), ylim=c(-4, 4), main=NA);
      text(-2, 3, paste('R = ', signif(cor(data[,16], data[, 5], use="na.or.complete"), d=2), sep=''), col='red')
      aabline(0, 1, lwd=1.,lty=2, col='red')
      #jj = which(!is.na(data[,16]) & !is.na(data[,8]))
      #fit = lm(data[jj,16]~data[jj,8]-1)
      #cd = cooks.distance(fit)
      #dev.off()
      
      pdf.name = paste("myplots/FIGURES/Influence_replicate_ZT21.pdf", sep='')
      pdf(pdf.name, width=6, height=3)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      
      cex = 0.1
      par(mfrow=c(1, 2))
      plot(data[which(rr<=cutoff), c(8, 16)], cex=cex, xlab='ZT21-Rep1', ylab='ZT21-Rep2', xlim=c(-4, 4), ylim=c(-4, 4), main=NA);
      points(data[which(rr>cutoff), c(8, 16)], cex=cex, col='orange');
      text(-2, 3, paste('R = ', signif(cor(data[which(rr<=cutoff) ,8], data[which(rr<=cutoff), 16], use="na.or.complete"), d=2), sep=''), col='red')
      abline(0, 1, lwd=1.5, col='red')
      
      jj = which(test$qval<0.05)
      plot(((test$phase.init[jj]-test$phase.refit[jj])%%24-24), test$pval.refit[jj], log='y', col='orange', cex=0.8, pch=1, 
           xlim = c(-12, 12), ylim =  c(10^-7, 1), xlab='phase diff with/without ZT21', ylab='pval without ZT21', axes=FALSE)
      axis(1,at=seq(-12, 12, by=6),cex.axis =1.0)
      axis(2, las=3,cex.axis = 0.8)
      abline(v=0, lwd=1.5, col='red', lty=2)
      abline(h=0.05, lwd=1.5, col='red', lty=2.0)
      box()
      #plot(test$phase.init[jj], test$phase.refit[jj], cex=1.0)
      #abline(0, 1, lwd=1.5, col='red')
      #abline(-24, 1, lwd=1.5, col='red')
      
      dev.off()
      
      #write.table(test, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Revision/Table_S0_outliers_ZT21.txt', 
      #            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
      #data = 2^data
      #pairs(data)
      #test = data.frame(xx$Gene.names[kk], rr, kk, stringsAsFactors = FALSE)
      #test = test[order(-test[,2]), ]
      
      #pdf.name = paste("myplots/FIGURES/Correlation_replicates_rhythmic_ones.pdf", sep='')
      #pdf(pdf.name, width=12, height=6)
      
      
      #dev.off()
      #jj = which(!is.na(data[,n]) & !is.na(data[, (n+8)]))
      #signif(cor(data[jj,n], data[jj, (n+8)]), d=2)
      
      #### check the rhythmicity and phases without ZT21
      data2 = as.matrix(xx[, c(5:20)])
      #rownames(data) = xx$Gene.names[kk]
      data2[,8] = NA;
      res = t(apply(data2, 1, f24_R2_alt2, t=c(0:15)*3))
      res[, 4] = t(apply(2^data2, 1, f24_R2_alt2, t=c(0:15)*3))[,4]
      res = cbind(res, qv=qvals(res[,6]))
      
      length(which(res[,7]<0.05))
      length(intersect(which(res[,7]<0.05), table.sx[,1]))
      length(setdiff(which(res[,7]<0.05), table.sx[,1]))
      length(setdiff(table.sx[,1], which(res[,7]<0.05)))
      
      index =table.sx[,1]
      
      pdf.name = paste("myplots/FIGURES/Influence_replicate_ZT21_to_rhythmic_proteins_1.pdf", sep='')
      pdf(pdf.name, width=1., height=1.)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      #par(mfrow=c(2, 2))
      cex = 0.1
      par(mgp = c(0.1,0.2,0.), mar = c(1.05,0.2,0.8,0.), pty='s');
      cex.main = 0.7;
      cex.axis = 0.6;
      xlab=NA;ylab=NA;
      lwd = 0.8
      plot(xx$pval.WT[index], res[index,6], cex=cex, log='xy', xlab=xlab, ylab=ylab, main='rhythmicity (p-val)', cex.main=cex.main,
           xlim=c(10^-10, 1), ylim=c(10^-10, 1), axes = FALSE)
      axis(1, at=c(10^-10, 10^-6, 10^-2), cex.axis = cex.axis)
      axis(2, at=c(10^-10, 10^-6, 10^-2), las=3,cex.axis = cex.axis)
      box()
      abline(0, 1, lwd=lwd, col='red')
      abline(h=0.05,lwd=lwd, col='red', lty=2)
      dev.off()
      #abline(h=0.01,lwd=1.5, col='red')
      
      pdf.name = paste("myplots/FIGURES/Influence_replicate_ZT21_to_rhythmic_proteins_2.pdf", sep='')
      pdf(pdf.name, width=1., height=1.)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      #par(mfrow=c(2, 2))
      cex = 0.1
      par(mgp = c(0.1,0.2,0.), mar = c(1.05,0.2,0.8,0.), pty='s');
      
      plot(xx$phase.WT[index], res[index,5], cex=cex, log='',xlab=xlab, ylab=ylab, main='phase', cex.main=cex.main, axes = FALSE)
      axis(1, at = seq(0, 24, by=6), cex.axis =cex.axis)
      axis(2, at = seq(0, 24, by=6), las=3,cex.axis = cex.axis)
      box()
      abline(0, 1, lwd=lwd, col='red')
      abline(-24, 1, lwd=lwd, col='red')
      abline(24, 1, lwd=lwd, col='red')
      dev.off()
      
      pdf.name = paste("myplots/FIGURES/Influence_replicate_ZT21_to_rhythmic_proteins_3.pdf", sep='')
      pdf(pdf.name, width=1., height=1.)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3, pty='s')
      #par(mfrow=c(2, 2))
      cex = 0.1
      par(mgp = c(0.1,0.2,0.), mar = c(1.05,0.2,0.8,0.), pty='s');
      plot(xx$relamp.WT[index], res[index,4], cex=cex, log='xy', xlab=xlab, ylab=ylab, main='rel.amp',cex.main = cex.main,
           xlim=c(0.01, 2), ylim=c(0.01, 2), axes = FALSE)
      abline(0, 1, lwd=lwd, col='red')
      axis(1, at = c(0.01, 0.1, 1), las=1, cex.axis =cex.axis)
      axis(2, at = c(0.01, 0.1, 1), las=3,cex.axis = cex.axis)
      box()
      dev.off()
      #abline(-24, 1, lwd=1.5, col='red')
      #abline(24, 1, lwd=1.5, col='red')
     
      
      #### corelation of mRNA across time points
      data = as.matrix(mrna[, c(2:13)])
      
      pdf.name = paste("myplots/FIGURES/Correlation_replicates_all_timepoints_mRNAs.pdf", sep='')
      pdf(pdf.name, width=6, height=4)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      
      par(mfrow=c(2, 3))
      for(n in 1:6)
      {
        #print(c((n-1)*3, (n+7)*3))
        plot(data[, c(n, (n+6))], cex=0.2, xlab=paste('ZT', ((n-1)*4+2), '-Rep1', sep = ''),
             ylab=paste('ZT', ((n-1)*4+2),  '-Rep2', sep=''), xlim=c(-4.5, 15), ylim=c(-4.5, 15), main=NA);
        text(-1,12, paste('R = ', signif(cor(data[,n], data[, (n+6)], use="na.or.complete"), d=3), sep=''), col='red')
        abline(0, 1, lwd=1.5, col='red')
      }

      dev.off()
      
    }
    #kk = table.sx[,1]
    Check.all.protein.4conditions = FALSE
    if(Check.all.protein.4conditions)
    {
        kk = grep('ZT', colnames(nuclear))
        boxplot(nuclear[,kk], cex=0.7, ylim=c(-2.5, 2.5), col=c(rep('white', 16), rep('red', 4), rep('gray', 4), rep('orange', 4)))
    }
    
}

#################
####### FIGURE 1 (Main + SUpplementary)
#################
FIGURE_1 = TRUE
if(FIGURE_1)
{
    #### Distribution of sample numbers
    pdf.name = paste("myplots/FIGURES/Distribution_time_points.pdf", sep='')
    pdf(pdf.name, width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    hist(nuclear$nb.timepoints, breaks=c(-1:16), col='gray', xlab=NA, ylab=NA, main=NA, axes=FALSE)
    axis(1,at=seq(0, 16, by=4),cex.axis =1.0)
    axis(2, at= seq(0, 3000, by=1000), las=1,cex.axis = 1.0)
    
    dev.off()
    
    #### localization comparison between nuclear and total proteins
    #### Use the GO Term annotation
    source('functions_nuclear.R')
    total = read.table('Tables_DATA/Table_total_proteins_PNAS.txt',sep='\t', header=TRUE)
    total.names = read.table('Tables_DATA/Table_total_proteins_index_names.txt', sep='\t', header=TRUE, as.is=c(2))
    
    cyto = read.delim('Annotations/uniprot-locations-Cytoplasm-mouse.tab', sep='\t', header=TRUE)
    nucleus = read.delim('Annotations/uniprot-locations-Nucleus-mouse.tab', sep='\t', header=TRUE)
    skeleton = read.delim('Annotations/uniprot-locations-Cytoskeleton-mouse.tab', sep='\t', header=TRUE)
    cyto = cyto$Gene.names
    cyto = unlist(strsplit(as.character(cyto), ' '))
    cyto = unique(cyto)
    nucleus = nucleus$Gene.names
    nucleus = unlist(strsplit(as.character(nucleus), ' '))
    nucleus = unique(nucleus)
    
    mix = intersect(nucleus, cyto)
    cyto = setdiff(cyto, mix)
    nucleus = setdiff(nucleus, mix)
    
    skeleton = unique(unlist(strsplit(as.character(skeleton$Gene.names), ' ')))
    length(intersect(skeleton, cyto))
    length(intersect(skeleton, nucleus))
    length(intersect(skeleton, mix))
    cytoskeleton = intersect(skeleton, cyto)
    
    
    nn = which(!is.na(match(nuclear.names[,3], nucleus)))
    mm = which(!is.na(match(nuclear.names[,3], mix)))
    cc = which(!is.na(match(nuclear.names[,3], cyto)))
    ss = which(!is.na(match(nuclear.names[,3], cytoskeleton)))
    
    index.nucleus = unique(nuclear.names[nn,1])
    index.mix = unique(nuclear.names[mm,1])
    index.cyto = unique(nuclear.names[cc,1])
    index.cytos = unique(nuclear.names[ss,1])
    index = unique(c(index.cyto, index.nucleus, index.mix, index.cytos))
    
    length(index)
    length(index.nucleus)
    length(index.mix)
    length(index.cyto)
    length(index.cytos)
    length(index.nucleus)/length(index)
    length(index.mix)/length(index)
    length(index.cyto)/length(index)
    length(index.cytos)/length(index)
    
    locas = rep(NA, nrow(nuclear))
    locas[index.cyto] = 'Cytoplasm'
    locas[index.cytos] = 'Cytoskeleton'
    locas[index.mix] = 'Nucleus/Cytoplasm'
    locas[index.nucleus] = 'Nucleus'
   
    xx = data.frame(nuclear.all, localization=locas, stringsAsFactors=FALSE)
    
    kk = grep('.total', colnames(xx))
    yy = xx[, c(1:45, kk, ncol(xx))]
    
    jj = which(!is.na(yy$gene.total)==TRUE)
    yy = yy[jj,]
    
    #write.table(yy, file='Tables_DATA/Nuclear_Total_prots_4Daniel.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    #x=x/s
    pdf.name = paste("myplots/FIGURES/localization_nuclear_prots.pdf", sep='')
    pdf(pdf.name, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    index = which(!is.na(xx$localization)==TRUE)
    index.nucleus = which(xx$localization=='Nucleus')
    index.mix = which(xx$localization=='Nucleus/Cytoplasm')
    index.cyto = which(xx$localization=='Cytoplasm')
    index.cytoskeleton = which(xx$localization=='Cytoskeleton')
    
    pie.sales= c(length(index.nucleus)/length(index), length(index.mix)/length(index), (length(index.cyto))/length(index), 
                 length(index.cytoskeleton)/length(index))
    names(pie.sales) <- c("47.4%", "27.3%","18.8%", "6.5%")
    cols = c("red","darkorange","deepskyblue","dodgerblue3")
    pie(pie.sales,col=cols, cex=1.1)
    #legend(0.65, 1.25, legend = c(' ',' ', ' ', ' ') , cex=0.9, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()
    
    #library(plotrix)
    #slices <- c(length(index.nucleus)/length(index), length(index.mix)/length(index), (length(index.cyto)-length(index.cytos))/length(index), length(index.cytos)/length(index))
    #lbls <- c("Nucleus", "Nucleus/Cytoplasm", "Cytoplasm", "Cytosekelton")
    #pie3D(slices,labels=lbls,explode=0.1, main="Pie Chart of Countries ")
    
    nn = which(!is.na(match(total.names[,2], nucleus)))
    mm = which(!is.na(match(total.names[,2], mix)))
    cc = which(!is.na(match(total.names[,2], cyto)))
    ss = which(!is.na(match(total.names[,2], cytoskeleton)))
    
    index.nucleus = unique(total.names[nn,1])
    index.mix = unique(total.names[mm,1])
    index.cyto = unique(total.names[cc,1])
    index.cytos = unique(total.names[ss,1])
    index = unique(c(index.cyto, index.nucleus, index.mix, index.cytos))
    
    length(index)
    length(index.nucleus)/length(index)
    length(index.mix)/length(index)
    length(index.cyto)/length(index)
    length(index.cytos)/length(index)
    
    locas = rep(NA, nrow(total))
    locas[index.cyto] = 'Cytoplasm'
    locas[index.cytos] = 'Cytoskeleton'
    locas[index.mix] = 'Nucleus/Cytoplasm'
    locas[index.nucleus] = 'Nucleus'
    
    xx = data.frame(total, localization=locas, stringsAsFactors=FALSE)
    #write.table(xx, file='Tables_DATA/Total_prots_4Daniel.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    pdf.name = paste("myplots/FIGURES/localization_total_prots.pdf", sep='')
    pdf(pdf.name, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    
    index = which(!is.na(xx$localization)==TRUE)
    index.nucleus = which(xx$localization=='Nucleus')
    index.mix = which(xx$localization=='Nucleus/Cytoplasm')
    index.cytoskeleton = which(xx$localization=='Cytoskeleton')
    index.cyto = which(xx$localization=='Cytoplasm')
    
    
    pie.sales= c(length(index.nucleus)/length(index), length(index.mix)/length(index), (length(index.cyto))/length(index), 
                 length(index.cytoskeleton)/length(index))
    
    names(pie.sales) <- c("25.3%", "30.2%","33.8%", "10.7%")
    cols = c("red","darkorange","deepskyblue","dodgerblue3")
    pie(pie.sales,col=cols, cex=1.1)
    #legend('topright', legend = c('','', '', '') , cex=1.0, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()
    
    ##########
    #### abosute intensities detected by MS for nuclear and total proteins
    pdf.name = paste("myplots/FIGURES/Intensity_localization_nuclear_prots.pdf", sep='')
    pdf(pdf.name, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    intense = c(160, 132, 18.2, 3.94)*10^10
    intense = intense/sum(intense)
    pie.sales= intense
    names(pie.sales) <- c("50.9%", "42.0%","5.8%", "1.3%")
    cols = c("red","darkorange","deepskyblue","dodgerblue3")
    pie(pie.sales,col=cols, cex=1.1)
    #legend('topright', legend = c('','', '', '') , cex=1.0, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()

    pdf.name = paste("myplots/FIGURES/Intensity_localization_total_prots.pdf", sep='')
    pdf(pdf.name, width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    intense = c(282, 238, 784, 87.6)*10^10
    intense = intense/sum(intense)
    pie.sales= intense
    names(pie.sales) <- c("20.3%", "17.1%","56.3%", "6.3%")
    cols = c("red","darkorange","deepskyblue","dodgerblue3")
    pie(pie.sales,col=cols, cex=1.1)
    #legend('topright', legend = c('','', '', '') , cex=1.0, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()
    
    ###### quality control of nuclear proteins and comparison with Kislinger et al. nuclear proteome
    Comparison.Uhlen.Kislinger = FALSE
    if(Comparison.Uhlen.Kislinger)
    {
        load(file='Annotations/Uhlen_Kislinger_Compartments_Uniprot_predicted_localization.Rdata')
        #cs = nuclear.names[cc, 3]
        
        #### Uniprot, Compartment, Uhlen
        #### percentages of proteins being considered as nuclear proteins according to different database
        length(index.nuclear.uniprot)/nrow(nuclear)
        length(index.nuclear.compartments)/(nrow(nuclear)-length(which(nuclear$Gene.names=='')))
        
        mm = match(nuclear.names[,3], uhlen[,1])
        mm = which(!is.na(mm))
        length(unique(nuclear.names[mm, 1]))/nrow(nuclear)
        overlap = length(unique(nuclear.names[mm, 1]))
        
        mm = match(nuclear.names[,3], uhlen[which(uhlen[,9]>0),1])
        mm = which(!is.na(mm))
        length(unique(nuclear.names[mm, 1]))/overlap
        
        #load(file='Annotations/Uhlen_Kislinger_Compartments_Uniprot_predicted_localization.Rdata')
        Add.coverage = FALSE
        if(Add.coverage)
        {
          #### Compartment
          known = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_knowledge_full.tsv',sep='\t',header=FALSE)
          text = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_textmining_full.tsv',sep='\t',header=FALSE)
          predict = read.delim('Annotations/Localization_COMPARTMENTS/mouse_compartment_predictions_full.tsv',sep='\t',header=FALSE)
          
          known = rbind(known[, c(1:4)], text[, c(1:4)], predict[, c(1:4)])
          #ll = unique(known[,4])
          #ll = ll[order(ll)]
          gg = (unique(known[,2]))
          mm = match(nuclear.names[,3], gg)
          length(which(is.na(mm)==TRUE))
          length(unique(nuclear.names[which(is.na(mm)==TRUE), 1]))
          
          1-length(unique(nuclear.names[which(is.na(mm)==TRUE), 1]))/nrow(nuclear)
          
          ### Ulhen
          mm1 = match(nuclear.names[,3], uhlen[,1])
          #length(which(!is.na(mm1)==TRUE))
          length(unique(nuclear.names[which(!is.na(mm1)==TRUE), 1]))
          
          length(unique(nuclear.names[which(!is.na(mm1)==TRUE), 1]))/nrow(nuclear)
          
          mm2 = match(nuclear.names[,3], uhlen[which(uhlen$nuclear>0),1])
          #length(which(!is.na(mm2)==TRUE))
          length(unique(nuclear.names[which(!is.na(mm2)==TRUE), 1]))
          
          length(unique(nuclear.names[which(!is.na(mm2)==TRUE), 1]))/nrow(nuclear)
          
        }
        
        coverages = c(61.5, 69.5, 87.4) 
        nucleus = c(length(index.nuclear.uniprot)/nrow(nuclear), 0.4923,
                    length(index.nuclear.compartments)/(nrow(nuclear)-length(which(nuclear$Gene.names==''))))
        
        nucleus = nucleus*100
        counts <- data.frame(coverages, nucleus)
        counts = as.matrix(t(counts))
        #colnames(counts) = c('')
        
        pdf.name = paste("myplots/FIGURES/Percentages_cyto_nuclear_Prots.pdf", sep='')
        pdf(pdf.name, width=1.8, height=1.7)
        par(cex = 0.65, las = 1, mgp = c(1.6,0.5,0), mar = c(5,3,2,0.8)+0.1, tcl = -0.3)
        
        barplot(counts, main=NA, col=c('white', "darkgray"), ylim=c(0, 100), beside =TRUE, cex.lab=1.0)
        #legend = c(' ', ' '), args.legend=list(x=12, y=1.1, bty='n', cex=1.0), beside=TRUE, las=2, cex.lab=1.0)
        #abline(h=0.15, col='red', lwd=1.5)
        
        dev.off()
        
        
        ### overlap between Kislinger et al. and our database
        require('VennDiagram')
        cols = c("gray","blue")
        kk = which(!is.na(kislinger$gene)==TRUE)
        #jj = which(!is.na(nuclear.names[,3])==TRUE)
        jj = which(!is.na(nuclear$Gene.names)==TRUE & nuclear$Gene.names!='')
        
        y0 = c(1:nrow(nuclear))
        #y1 = c(1:nrow(kislinger))+5000
        mm = match(kislinger$gene, nuclear.names[,3])
        jj = which(!is.na(mm)==TRUE)
        ii = which(is.na(mm)==TRUE)
        #y1[jj] = nuclear.names[mm[jj], 1];
        y1 = c(c(1:length(jj)), c(1:length(ii))+5000)
        #venn.diagram(list(Kislinger = kislinger$gene[kk], ours= nuclear$Gene.names[jj,3]),
        venn.diagram(list(Kislinger = y1, ours= y0),
        height = 2.0, width = 2.0, units = "in",
        fill = cols, cat.col='black', resolution=500,
        alpha = c(0.6, 0.6), cex = 1.0, cat.fontface = 2,lty=1,sub.cex=0.8,
        mar = rep(0., 4)+0.1, cat.dist=c(0.15, 0.08), ext.dist=0.1,
        filename = "myplots/FIGURES/Kinslinger_nuclear_overlap.png");
    }
    
    #### Coverages according to RNA-seq and Uniprot
    load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/Uniprot_proteins_localizations.Rdata')
    length(intersect(nuclear.names[,3], gene.expressed))/length(gene.expressed)
    length(intersect(nuclear.names[,3], intersect(gene.expressed, nucleus)))/length(intersect(gene.expressed, nucleus))
    length(intersect(nuclear.names[,3], intersect(gene.expressed, mix)))/length(intersect(gene.expressed, mix))
    
    #length(which(!is.na(match(nuclear.names[,3], intersect(nucleus, gene.expressed)))))/length(intersect(nucleus, gene.expressed))
    #length(which(!is.na(match(nuclear.names[,3], mix))))/length(intersect(mix, gene.expressed))
    
    #### subcompartment of nuclear proteins
    subs = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Intra_nucleus_localisation_data.tsv', header=FALSE, sep='\t', as.is=c(1, 2))
    ids = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/uniprot-GeneIDs-mapping-human-mouse-rat.txt', header=TRUE, sep='\t', as.is=c(1:5))
    
    comparts = unique(unlist(strsplit(as.character(subs[,2]), ' ')))
    comparts = comparts[which(comparts!='')]
    
    kk = which(subs[,2]!='')
    subs = subs[kk,]
    
    mm = match(subs[,1], ids[,1])
    
    #poly = read.table('/Users/jiwang/RNA_seq_Data/Total-RNA-seq/mRNAs_pre_mRNAs_PolyA_RNA_seq_Gachon.txt', header=TRUE, sep='\t') ### this polyA RNA-seq data needs to be checked
    #ensgene = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/ensgenes_mapping.txt')
    #cutoff.expressed = 0
    #gene.expressed = unique(c(as.character(poly$gene[which(poly$mean.ex>=cutoff.expressed)]), as.character(nuclear.names[,3])))
    #gene.expressed = unique(unlist(strsplit(as.character(gene.expressed), ';')))
    
    refs = c()
    for(n in 1:nrow(subs))
    {
        cat(n, '\n')
        prot = as.character(subs[n, 1])
        mm = match(prot, as.character(ids[,1]))
        mm = mm[which(!is.na(mm)==TRUE)]
        
        locas = subs[n, 2]
        locas = unlist(strsplit(as.character(locas), ' '))
        nb.locas = length(locas)
        
        if(length(mm)==0){### no gene symbol found
            for(loca in locas)
            {
                refs = rbind(refs, c(prot, loca, NA, NA, NA))
            }
        }else{ ### found gene symbols
            
            symbols = ids[mm, 5]
            symbols = unlist(strsplit(as.character(symbols), ' '))
            
            jj = match(symbols, gene.expressed)
            jj = jj[which(!is.na(jj)==TRUE)]
            #kk = match(symbols, )
            kk = match(symbols, nuclear.names[,3])
            kk = kk[which(!is.na(kk)==TRUE)]
            #kk = which(!is.na(jj)==TRUE)
            
            for(loca in locas)
            {
                refs = rbind(refs, c(prot, loca, symbols[1], length(jj)>0, length(kk)>0))
            }
        }
    }
    
    refs = data.frame(refs, stringsAsFactors=FALSE)
    
    colnames(refs) = c('protein.ID', 'location', 'gene.symbol', 'expressed', 'detected.nuclear.prot')
    
    comparts = unique(refs[,2])
    comparts = comparts[which(comparts!='')]
    counts = c()
    for(compart in comparts)
    {
        #print(length(which(refs[,2]==compart)))
        #print(length(which(refs[,2]==compart & refs[,4]==1)))
        #print(length(which(refs[,2]==compart & refs[,4]==1))/length(which(refs[,2]==compart)))
        #print(length(which(refs[,2]==compart & refs[,4]==1 & refs[,5]==1))/length(which(refs[,2]==compart & refs[,5]==1)))
        counts = rbind(counts, c(length(which(refs[,2]==compart)), length(which(refs[,2]==compart & refs[,4]==TRUE)),
        length(which(refs[,2]==compart & refs[,4]==TRUE & refs[,5]==TRUE)),
        length(which(refs[,2]==compart & refs[,4]==TRUE & refs[,5]==TRUE))/length(which(refs[,2]==compart & refs[,4]==TRUE))))
    }
    
    counts = data.frame(comparts, counts, stringsAsFactors=FALSE)
    colnames(counts) = c('compartments', 'nb.prots','nb.prots.expressed', 'prots.detected.nuclear', 'percentages')
    
    pdf.name = paste("myplots/FIGURES/Subcompartments_Nuclear_Prots.pdf", sep='')
    pdf(pdf.name, width=2.5, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(5,3,2,0.8)+0.1, tcl = -0.3)
    
    percents = c(length(intersect(nuclear.names[,3], gene.expressed))/length(gene.expressed), length(intersect(nuclear.names[,3], intersect(gene.expressed, nucleus)))/length(intersect(gene.expressed, nucleus)), counts[,5])*100
    barplot(percents, main=NA, col=c(rep("white", 1), rep("bisque1",1), rep("azure1", 8)), ylim=c(0, 100), space=c(0.1, 0.6, 0.6, rep(0.2, 7)),
    beside=TRUE, las=2, cex.lab=1.0)
    #abline(h=1.0, col='black', lwd=2.0)
    dev.off()
    
    #############
    #### Coverages and rhythmicity percentages for different functions
    ############
    functional.annotations = FALSE
    if(functional.annotations)
    {
        #xx = read.delim('Annotations/kinases_phosphatase/uniprot-kinase.tab',header=TRUE, sep='\t', as.is = c(3));
        ####### Annotations from GOTERMS TFs, coregulators, RNA-binding proteins, kinase, phosphatase
        mapping.all = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Manually_curation_TFs_pools.txt', header=TRUE, sep='\t')
        tfs = mapping.all$TFs.pool
        tfs = unlist(strsplit(as.character(tfs), ';'))
        mcors = c('Ncor1', 'Ep300', 'Hdac1', 'Hdac2', 'Hdac5', 'Sirt6', 'Crebbp', 'Ncoa1', 'Ncoa2', 'Ncoa3', 'Ncor1', 'Ncor2', 'Nrip1', 'Kat2a', 'Kat2b', 'Hdac1', 'Hdac2',
        'Hdac3', 'Hdac4', 'Hdac5', 'Hdac6', 'Hdac7', 'Hdac8', 'Hdac9', 'Crtc1', 'Crtc2', 'Crtc3', 'Rb1', 'Med1', 'Ppargc1a', 'Ppargc1b', 'Setd2', 'Setdb1')
        tfs = tfs[which(is.na(match(tfs, mcors))==TRUE)]
        
        cofactor = read.delim('Annotations/Cofactors_goterms.txt',header=FALSE, sep='\t', as.is = c(3));
        cofactor = unique(cofactor[,3])
        
        rna.processing = rna.processing.protein
        kinase = read.delim('Annotations/kinases_phosphatase/kinase_goterm.txt',header=FALSE,sep='\t', as.is = c(3))
        kinase = unique(kinase[,3])
        phosphatase = read.delim('Annotations/kinases_phosphatase/phosphatase_goterm.txt',header=FALSE, sep='\t', as.is=c(3));
        phosphatase = unique(phosphatase[,3])
        
        dna.repair = read.delim('Annotations/DNA_repair_goterm.txt', header=FALSE, sep='\t', as.is=c(3));
        dna.repair = unique(dna.repair[,3])
        
        rRNA = read.delim('Annotations/ribosome_biogenesis_goterm.txt', header=FALSE, sep='\t', as.is=c(3));
        rRNA = unique(rRNA[,3])
        
        load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/Uniprot_proteins_localizations.Rdata')
        load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/List_genes_expresed_mouse_liver.Rdata')
        #poly = read.table('/Users/jiwang/RNA_seq_Data/Total-RNA-seq/mRNAs_pre_mRNAs_PolyA_RNA_seq_Gachon.txt', header=TRUE, sep='\t') ### this polyA RNA-seq data needs to be checked
        #ensgene = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/ensgenes_mapping.txt')
        #cutoff.expressed = 0
        #gene.expressed = unique(c(as.character(poly$gene[which(poly$mean.ex>=cutoff.expressed)]), as.character(nuclear.names[,3])))
        #gene.expressed = unique(unlist(strsplit(as.character(gene.expressed), ';')))
        
        ### Localization annotations
        #cyto = read.delim('Tables_DATA/Annotation_pure_cytoplasmic_proteins_v2.txt', sep='\t', header=TRUE, as.is=TRUE)
        #nucleus = read.delim('Tables_DATA/Annotation_pure_nucleus_proteins_v2.txt', sep='\t', header=TRUE, as.is=TRUE)
        #mix = read.delim('Tables_DATA/Annotation_mix_cytoplasmic_nucleus_proteins_v2.txt', sep='\t', header=TRUE, as.is=TRUE)
        
        tfs.expressed = tfs[which(!is.na(match(tfs, gene.expressed))==TRUE)]
        cors.expressed = cofactor[which(!is.na(match(cofactor, gene.expressed))==TRUE)]
        rna.processing.expressed = rna.processing[which(!is.na(match(rna.processing, gene.expressed))==TRUE)]
        
        
        kinase.expressed = kinase[which(!is.na(match(kinase, gene.expressed))==TRUE)]
        kinase.localization = rep(NA, length(kinase.expressed))
        kk = match(kinase.expressed, cyto)
        kinase.localization[which(!is.na(kk))] = 'Cytoplasm'
        kk = match(kinase.expressed, mix)
        kinase.localization[which(!is.na(kk))] = 'Nucleus/Cytoplasm'
        kk = match(kinase.expressed, nucleus)
        kinase.localization[which(!is.na(kk))] = 'Nucleus'
        kk = match(kinase.expressed, nuclear.names[,3])
        index = nuclear.names[kk, 1]
        
        keep.kinase = data.frame(kinase.expressed, kinase.localization, index, stringsAsFactors=FALSE)
        jj = which(keep.kinase[,2]=='Cytoplasm')
        ii = which(keep.kinase[,2]!='Cytoplasm' & !is.na(keep.kinase[,3])==TRUE)
        
        phosphatase.expressed = phosphatase[which(!is.na(match(phosphatase, gene.expressed))==TRUE)]
        phosphatase.localization = rep(NA, length(phosphatase.expressed))
        kk = match(phosphatase.expressed, cyto)
        phosphatase.localization[which(!is.na(kk))] = 'Cytoplasm'
        kk = match(phosphatase.expressed, mix)
        phosphatase.localization[which(!is.na(kk))] = 'Nucleus/Cytoplasm'
        kk = match(phosphatase.expressed, nucleus)
        phosphatase.localization[which(!is.na(kk))] = 'Nucleus'
        kk = match(phosphatase.expressed, nuclear.names[,3])
        index = nuclear.names[kk, 1]
        
        keep.phosphatase = data.frame(phosphatase.expressed, phosphatase.localization, index, stringsAsFactors=FALSE)
        
        #length(kinase.expressed)
        #nrow(keep.kinase)
        #length(which(!is.na(keep.kinase[,3])==TRUE))
        #length(which(keep.kinase[,2]=='Cytoplasm'))
        #length(which(keep.kinase[,2]=='Cytoplasm' & !is.na(keep.kinase[,3])==TRUE))
        length(which(keep.kinase[,2]=='Nucleus'))
        length(which(keep.kinase[,2]=='Nucleus' & !is.na(keep.kinase[,3])==TRUE))
        length(which(keep.kinase[,2]=='Nucleus/Cytoplasm'))
        length(which(keep.kinase[,2]=='Nucleus/Cytoplasm' & !is.na(keep.kinase[,3])==TRUE))


        length(which(keep.phosphatase[,2]=='Nucleus'))
        length(which(keep.phosphatase[,2]=='Nucleus' & !is.na(keep.phosphatase[,3])==TRUE))
        length(which(keep.phosphatase[,2]=='Nucleus/Cytoplasm'))
        length(which(keep.phosphatase[,2]=='Nucleus/Cytoplasm' & !is.na(keep.phosphatase[,3])==TRUE))
        
        #save(keep.kinase, keep.phosphatase, file='Rdata/kinases_phosphatases_keep.Rdata')
        ##### Statistics for different functions
        stat = c()
        
        xx = tfs.expressed
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        xx = cors.expressed
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        xx = rna.processing.expressed
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        xx = kinase.expressed
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        xx = intersect(kinase.expressed, nucleus)
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        
        xx = phosphatase.expressed
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))
        
        xx = intersect(phosphatase.expressed, nucleus)
        kk = match(xx, nuclear.names[,3])
        kk = kk[which(!is.na(kk))]
        index = nuclear.names[kk, 1]
        index = index[which(nuclear$qv[index]<0.05)]
        stat = rbind(stat, c(length(kk)/length(xx), length(index)/length(kk)))
        print(c(length(kk)/length(xx), length(index)/length(kk)))

        
        stat = data.frame(c('TFs', 'CoRs','RNA.Processing', 'Kinases', 'Kinases', 'Phosphatases',  'Phosphatases'), stat, stringsAsFactors=FALSE)
        colnames(stat) = c('functions', 'coverage', 'rhythmic.percent')
        save(stat, file='Rdata/functionality_coverage_rhythmic_percentages_nuclear_prots.Rdata')
        
    }
    
    load(file='Rdata/functionality_coverage_rhythmic_percentages_nuclear_prots.Rdata')
    type= c('TFs', 'CoRs','RNA.Processing', 'Kinases', 'Kinases','Phosphatases' ,'Phosphatases')
    #coverages = stat[,2]*(stat[,3])
    #rhythm = stat[,2]*(1-stat[,3])
    #counts <- rbind(coverages, rhythm)*100
    counts = stat[,2]*100
    
    pdf.name = paste("myplots/FIGURES/Coverage_Nuclear_Prots.pdf", sep='')
    pdf(pdf.name, width=2.5, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(5,3,2,0.8)+0.1, tcl = -0.3)
    barplot(counts, main=NA, col=c(rep("gray", 4), 'gray40', 'gray', 'gray40'), ylim=c(0, 100), space=c(0.1, rep(0.5, 3), 0.15, 0.5, 0.15),las=2, cex.lab=1.0,
    legend = c(' ', ' '), args.legend=list(fill=c('gray', 'gray40'), x=7, y=120, bty='n', cex=1.0))
    #abline(h=0.15, col='red', lwd=1.5)
    dev.off()
    
    ########
    ##### Summary of clock genes detected in nuclear data
    #########
    examples = c('Arntl', 'Clock', 'Per1', 'Per2','Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3')
    mm = match(examples, nuclear$Gene.names)
    phase = nuclear$phase[mm]
    fch = nuclear$amp[mm]
    cbind(examples, phase, fch)
    
    names = examples;
    names[1] = ''
    names[2] = 'Arntl,Clock'
    norm = fch;
    #phase=phase_comp
    #names=	names_comp
    o1 = order(phase)
    phase = phase[o1]
    names = names[o1]
    norm = norm[o1]
    
    ### modify the phase gaps in order to better represent phases
    for(n in 1:(length(phase)-1))
    {
        diff = abs(phase[n+1] - phase[n])
        if(diff<0) diff = diff+24
        if(diff>=24) diff = diff-24
        if(diff<0.5) phase[n+1] = phase[n]+0.5
        if(phase[n+1]>=24) phase[n+1] = phase[n+1]-24
    }
    
    rr = ceiling(max(norm))
    rr = c(-rr, rr)
    ampl1 = norm
    aa = ampl1*cos(2*pi/24*phase)
    bb = ampl1*sin(2*pi/24*phase)
    CC= (aa -1i*bb) * exp(1i*pi/2)
    a = Re(CC)
    b = Im(CC)
    
    ampl1 = ampl1+0.3
    
    aa = ampl1*cos(2*pi/24*phase)
    bb = ampl1*sin(2*pi/24*phase)
    
    CC= (aa -1i*bb) * exp(1i*pi/2)
    aa = Re(CC)
    bb = Im(CC)
    
    ### tune the distance between clock and bmal1
    #phase[2] = 4.5
    #phase[3] = 6.5
    rr = c(-5,5)
    
    pdf.name = paste("myplots/FIGURES/Summary_Clock_Genes.pdf", sep='')
    pdf(pdf.name, width=3.8, height=3.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    plot(a, b, main=NA, type='n', xlim=rr, ylim=rr, axes=FALSE, xlab='', ylab='', lwd=1.6, pch=21, col='darkblue',bg='red',cex=1.0,bty='n')
    #abline(v=0,h=0, col='darkgray',lwd=2.0)
    arrows(-5, 0, 5, 0,  col='darkgray', lty=2, lwd=1., length=0.0);
    arrows(0, -5, 0, 5, col='darkgray', lty=2, lwd=1., length=0.0);
    
    phi=seq(0,2*pi,len=1000)
    lwd = 1.0
    lines(5*cos(phi), 5*sin(phi), col='darkgray', lwd=2.0)
    lines(4*cos(phi), 4*sin(phi), col='darkgray', lty=2,lwd=lwd)
    lines(3*cos(phi), 3*sin(phi), col='darkgray', lty=2, lwd=lwd)
    lines(2*cos(phi), 2*sin(phi), col='darkgray', lty=2,lwd=lwd)
    lines(1*cos(phi), 1*sin(phi), col='darkgray', lty=2,lwd=lwd)
    
    for(n in 1:length(phase))
    {
        col='black'
        bg.col = 'red'
        if(names[n]=='CK1'|names[n]=='AMPK'|names[n]=='GSK3') bg.col = 'orange';
        if(names[n]=='Nr1d1'|names[n]=='Nr1d2'|names[n]=='Rora'|names[n]=='Rorc') bg.col = 'green';
        if(names[n]=='Dbp'|names[n]=='Tef'|names[n]=='Hlf'|names[n]=='Nfil3'|names[n]=='Bhl40') bg.col = 'orange';
        points(a[n],b[n], type='p',lwd=1.6, pch=21, col=col,bg=bg.col, cex=1.2,bty='n')
        
        
        pos = 3;
        if(names[n]=='Arntl,Clock') pos=3;
        if(names[n]=='Per1') pos=2;
        if(names[n]=='Per2') pos=2;
        if(names[n]=='Per3') pos=2;
        if(names[n]=='Rorc') pos=1;
        if(names[n]=='Rora') pos=4;
        if(names[n]=='CK1') pos=2;
        if(names[n]=='AMPK') pos=2;
        if(names[n]=='GSK3') pos=4;
        if(names[n]=='Nfil3') pos=1;
        if(names[n]=='Cry2') pos=1;
        
        text(a[n], b[n], toupper(names[n]), cex=0.8,col=col, pos=pos, offset=0.5)
    }
    legend('topright', legend = c(' ',' ', ' ') , cex=1.2, pch=c(21, 21,21), col='black', pt.cex=1.2, pt.lwd=1.0, pt.bg= c('red', 'green', 'orange'), border = NA, bty = 'n')
    
    dev.off()
    
    
    ###########
    ### Individal expamples of core clock genes and clock-related genes
    ###########
    require('plotrix')
    examples = c('Arntl', 'Clock', 'Per1', 'Per2', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp',  'Tef', 'Hlf', 'Nfil3')
    mm = match(examples, nuclear$Gene.names)
    for(gene in examples)
    {
        pdf.name = paste("myplots/FIGURES/Fig_examples_WT_KO_not_centered_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.5, height=1.5)
        par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        n = which(nuclear$Gene.names==gene)
        if(length(n)>1) n = n[which(n==min(n))]
        name = nuclear[n,20]
        if(name=='Csnk1d;Csnk1e') name = 'Csnk1d/e'
        if(name=='Prkab1;Prkab2') name = 'Prkab1/2'
        y01 = as.numeric(nuclear[n,c(1:8)])
        y02 = as.numeric(nuclear[n,c(9:16)])
        y0 = rbind(y01, y02)
        
        averg = c()
        err = c()
        for(ii in 1:8)
        {
            xx = y0[,ii]
            kk = length(which(is.na(xx)==TRUE))
            if(kk==0){averg = c(averg, mean(xx)); err=c(err, sd(xx)/sqrt(2))}
            if(kk==1){averg = c(averg, xx[which(!is.na(xx)==TRUE)]); err=c(err, 0)}
            if(kk==2){averg = c(averg, NA); err=c(err, NA)}
            
        }
        #y0 = (y0-mean(y0[which(!is.na(y0)==TRUE)]))
        kk = grep('Cry.KO', colnames(nuclear))
        y1 = as.numeric(nuclear[n, kk])
        #if(length(which(!is.na(y1))==TRUE)>0) y1 = (y1 - mean(y1[which(!is.na(y1)==TRUE)]))
        
        #kk = grep('Bmal.WT', colnames(nuclear))
        #y2 = as.numeric(nuclear[n, kk])
        #if(length(which(!is.na(y2))==TRUE)>0) y2 = (y2 - mean(y2[which(!is.na(y2)==TRUE)]))
        #kk = grep('Bmal.KO', colnames(nuclear))
        #y3 = as.numeric(nuclear[n, kk])
        #if(length(which(!is.na(y3))==TRUE)>0) y3 = (y3 - mean(y3[which(!is.na(y3)==TRUE)]))
        
        time = c(0:7)*3
        lims = range(c(averg-err, averg+err, y1), na.rm=TRUE)
        if(gene=='Cry1' | gene=='Cry2')
        lims = range(c(averg-err, averg+err), na.rm=TRUE)
        #lims = range(c(y0[which(!is.na(y0)==TRUE)]))
        
        plotCI(time, averg, err, scol="darkblue",lwd=1.2, pch=16, col="darkblue", main=toupper(name), ylim=lims, ylab=NA, cex=1.0, cex.lab=1.0, cex.main=1.2, xlab=NA, xlim=c(0,24), axes=FALSE)
        points(time, averg, col="darkblue", lwd=1.2, type='l')
        if(gene!='Cry1' & gene!='Cry2')
        points(6*c(0:3), y1, type='b', pch=16, lwd=1.2, cex=1.0, col="darkred")
       
        axis(1,at=6*c(0:4),cex.axis =1.2)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = signif(lims, d=1)
        by = signif((lims[2]-lims[1])/4,d=1)
        print(gene)
        print(lims)
        if(gene=='Arntl') axis(2,at = c(-0.6, -0.3, 0.0, 0.3),las=1,cex.axis = 1.2)
        if(gene=='Clock') axis(2,at = c(-0.5, -0.2, 0.1, 0.4),las=1,cex.axis = 1.2)
        if(gene=='Per1') axis(2,at = c(-2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Per2') axis(2,at = c(-2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Cry1') axis(2,at = c(-2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Cry2') axis(2,at = c(-1, -0.5, 0, 0.5),las=1,cex.axis = 1.2)
        if(gene=='Rora') axis(2,at = c(-2, -1.5, -1, -0.5, 0),las=1,cex.axis = 1.2)
        if(gene=='Rorc') axis(2,at = c(-3, -2, -1, -0, 1),las=1,cex.axis = 1.2)
        if(gene=='Nr1d1') axis(2,at = c(-3, -2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Nr1d2') axis(2,at = c(-2, -1, 0, 1, 2),las=1,cex.axis = 1.2)
        if(gene=='Dbp') axis(2,at = c(-4, -3, -2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Tef') axis(2,at = c(-3, -2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Hlf') axis(2,at = c(-4, -3, -2, -1, 0, 1),las=1,cex.axis = 1.2)
        if(gene=='Nfil3') axis(2,at = c(-1, 0, 1, 2),las=1,cex.axis = 1.2)
        
        box()
        
        abline(h=0,lty=2,lwd=1.5, col="darkgray")
        #legend(15,lims[2],c('WT','CRYDKO'),pch=c(16,16),lty=c(1,1),col=c("darkblue","darkred"))
        
        dev.off()
    }

}

#################
####### FIGURE 2 (Main + SUpplementary)
#################
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
#table.sx$nb.cryko = res.ko[,4]
#table.sx$cor.cryko = res.ko[,5]

#table.sx$pval.bmalwt = res.ko[,4]
#table.sx$nb.bmalwt = res.ko[,5]
#table.sx$cor.bmalwt = res.ko[,6]

#table.sx$pval.bmalko = res.ko[,7]
#table.sx$nb.bmalko = res.ko[,8]
#table.sx$cor.bmalko = res.ko[,9]


Table.Sup = FALSE
if(Table.Sup)
{
    index = table.sx[,1]
    xx = nuclear[index, ]
    xx = xx[, c(1:16, 34:37)]
    xx = data.frame(table.sx[, -c(1, 12, 13, 19)], xx, stringsAsFactors=FALSE)
    colnames(xx)[3:8] = paste(colnames(xx[3:8]), '.WT', sep='')
    #colnames(xx)[17:19] = c('chow.test.WT.KO', 'nb.timepoints.KO', 'cor.WT.KO')
    xx = xx[, c(1:2, 9:16, 3:8, 17:37)]
    #index = c(1:16, 34:37)
    #yy = table.sx
    #yy = yy[, c(1:3, 10:20, 4:9)]
    #xx = yy
    #colnames(xx)[15:20] = paste(colnames(xx)[15:20], '.WT', sep='')
    #xx = cbind(yy, nuclear[xx[,1], index])
    #xx = xx[, -1]
    write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S2_nuclear_proteins_Rhythmic.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

    
    index = table.sx[,1]
    xx = nuclear.all[index, ]
    
    kk = c(c(20, 19, 17, 18, 1:16, 34:37, 22:28), grep('mRNA', colnames(xx)), grep('mrna', colnames(xx)), grep('total', colnames(xx))[c(1:16, 18:24)])
    
    xx = xx[,kk]
    colnames(xx)[c(68, 74)] = c('nb.timepoints.total', 'qv.total')
    #colnames(xx)[25:31] = paste(colnames(xx)[25:31], '.WT', sep='')
    write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S2_nuclear_proteins_rhythmic_mRNA_Total.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    
}
FIGURE_2 = TRUE
if(FIGURE_2)
{
    ######
    ###Phases and amplitudes distributions for all rhythmic proteins detected
    ######
    qq = c(0:100)/100
    nb.rhythmic = c()
    for(n in 1:length(qq))
    {
        cutoff.qq = qq[n]
        nb.rhythmic = c(nb.rhythmic,length(which(nuclear$qv<cutoff.qq)))
    }
    
    pdf('myplots/FIGURES/Nb_rhythmic_Proteins_FDR.pdf', width=1.8, height=1.8)
    par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(qq, nb.rhythmic+0.01, type='l', lty=1, lwd=1.2, main=NA, log='y',ylim=c(1,5000), xlab=NA, ylab=NA, col=c('black'), axes=FALSE)
    #points(qq, nb.rhythmic.robles+0.01, type='l', lty=1, lwd=2, col='green')
    abline(v=0.05,col='darkblue',lwd=1.5,lty=2)
    axis(1,at=seq(0, 1.0, by=0.2),cex.axis =0.8)
    axis(2, at= c(1, 5, 50, 500, 5000), las=0,cex.axis = 0.8)
    box()
    #abline(v=0.2,col='darkgreen',lwd=2,lty=2)
    #abline(v=0.1,col='darkgreen',lwd=2,lty=2)
    #legend('topright', legend = c('Mauvoisin','Robles'), lty=c(1,1), cex=0.7,col = c('blue', 'green'), border = NA, bty = 'n')
    dev.off()
    
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
    
    #### All detected rhythmic proteins
    jj = which(nuclear$qv<0.05)
    phases = nuclear$phase[jj]
    amps = nuclear$amp[jj]
    
    pdf('myplots/FIGURES/Phase_distribution_all_rhythmic_detected.pdf', width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    pdf('myplots/FIGURES/Amplitudes_distribution_all_rhythmic_detected.pdf', width=1.7, height=1.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=c(0:10)/2;
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    lwd = 1.
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd)
    axis(1,at=seq(0, 5, by=1),cex.axis =1.0)
    axis(2, at= c(1, 5, 10, 50, 100, 200), las=1,cex.axis = 1.0)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd)
    lines(h$breaks, rep(0,length(h$breaks)), type='S')
    invisible(h)
    
    dev.off()
    
    #########
    ###Phases and amplitudes distributions for rhythmic proteins grouped with localizations and their comparisons with mRNAs
    ########
    Use.Compartment = FALSE
    if(!Use.Compartment)
    {
        library(plotrix)
        library("circular")
        make_circ_coord = function(t,x,ttot=24)
        {
            dt=(t[2]-t[1])*.45
            a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
            h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
            list(angles=a,heights=h)
        }
        
        #### Nuclear proteins
        jj = which(table.sx$Localization=='Nucleus')
        kk = which(table.sx$Localization=='Nucleus/Cytoplasm')
        
        phases.1 = table.sx$phase[jj]
        amps.1 = table.sx$amp[jj]
        phases.2 = table.sx$phase[kk]
        amps.2 = table.sx$amp[kk]
        
        pdf('myplots/FIGURES/Phase_Nucleus_both.pdf', width=2.2, height=2.2)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        #hist(phases,breaks=c(0:12)*2,xlab=NA,ylab=NA, main=NA,col='gray50',cex.lab=1.0,xlim=c(0,24),axes=FALSE)
        #axis(1,at=seq(0,24,by=4),cex.axis =0.7)
        #axis(2,at = seq(0,30, by=10),las=1,cex.axis = 0.7)
        cex.axis=0.7; cex.lab=0.1; lwd=0.5;
        
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(phases.1, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        col=rgb(0.2,1.0,0.2);
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=col)
        #circular_phase24H_histogram(phases, , )
        
        #cex.axis=0.75; cex.lab=0.1; lwd=0.7;
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(phases.2, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        col=rgb(0.1,0.6,0.1);
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=col, add=TRUE)
        
        
        dev.off()
        
        pdf('myplots/FIGURES/Figure_Amplitudes_Nuclear.pdf', width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        col=rgb(0.2,1.0,0.2);
        breaks=c(0:10)/2;
        h = hist(amps.1, breaks=breaks, plot=FALSE)
        y = h$counts
        y[which(y<=0)] = 0.5;
        
        lwd = 1.5
        plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd, col=col)
        axis(1)
        axis(2)
        lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd, col=col)
        lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd, col=col)
        lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd, col=col)
        lines(h$breaks, rep(0,length(h$breaks)), type='S', lwd=lwd, col=col)
        invisible(h)
        
        dev.off()
        
        pdf('myplots/FIGURES/Figure_Amplitudes_Both.pdf', width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        col=rgb(0.1,0.6,0.1);
        breaks=c(0:10)/2;
        h = hist(amps.2, breaks=breaks, plot=FALSE)
        y = h$counts
        y[which(y<=0)] = 0.5;
        lwd = 1.5;
        plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=1.0, col=col)
        #lines(h$breaks, c(NA,y), type='S', lwd=1.0, col=col)
        axis(1)
        axis(2)
        lines(h$breaks, c(h$counts,NA), type='s', lwd=1.0, col=col)
        lines(h$breaks, c(NA,h$counts), type='h', lwd=1.0, col=col)
        lines(h$breaks, c(h$counts,NA), type='h',lwd=1.0, col=col)
        lines(h$breaks, rep(0,length(h$breaks)), type='S', col=col)
        invisible(h)
        
        dev.off()
        
        #### Cytoplasmic and Cytoskeleton proteins
        jj = which(table.sx$Localization=='Cytoplasm')
        kk = which(table.sx$Localization=='Cyto.Cytoskeleton')
        
        phases.1 = table.sx$phase[jj]
        amps.1 = table.sx$amp[jj]
        phases.2 = table.sx$phase[kk]
        amps.2 = table.sx$amp[kk]
        
        pdf('myplots/FIGURES/Phase_Cyto_Cytoskeleton.pdf', width=2.2, height=2.2)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        #hist(phases,breaks=c(0:12)*2,xlab=NA,ylab=NA, main=NA,col='gray50',cex.lab=1.0,xlim=c(0,24),axes=FALSE)
        #axis(1,at=seq(0,24,by=4),cex.axis =0.7)
        #axis(2,at = seq(0,30, by=10),las=1,cex.axis = 0.7)
        cex.axis=0.7; cex.lab=0.1; lwd=0.7;
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(phases.1, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        col=rgb(1.0,0.4,0.5);
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=col)
        #circular_phase24H_histogram(phases, , )
        
        cex.axis=0.75; cex.lab=0.1; lwd=0.7;
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(phases.2, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        col=rgb(0.6,0,0.2)
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=col, add=TRUE)
        
        
        dev.off()
        
        pdf('myplots/FIGURES/Figure_Amplitudes_Cyto.pdf', width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        col=rgb(1.0,0.4,0.5);
        breaks=c(0:10)/2;
        h = hist(amps.1, breaks=breaks, plot=FALSE)
        y = h$counts
        y[which(y<=0)] = 0.5;
        
        lwd = 1.5
        plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=lwd, col=col)
        axis(1)
        axis(2)
        lines(h$breaks, c(h$counts,NA), type='s', lwd=lwd, col=col)
        lines(h$breaks, c(NA,h$counts), type='h', lwd=lwd, col=col)
        lines(h$breaks, c(h$counts,NA), type='h',lwd=lwd, col=col)
        lines(h$breaks, rep(0,length(h$breaks)), type='S', lwd=lwd, col=col)
        invisible(h)
        
        dev.off()
        
        pdf('myplots/FIGURES/Figure_Amplitudes_Cytoskeleton.pdf', width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        col=rgb(0.6,0,0.2)
        breaks=c(0:10)/2;
        h = hist(amps.2, breaks=breaks, plot=FALSE)
        y = h$counts
        y[which(y<=0)] = 0.5;
        lwd = 1.5;
        plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=1.0, col=col)
        #lines(h$breaks, c(NA,y), type='S', lwd=1.0, col=col)
        axis(1)
        axis(2)
        lines(h$breaks, c(h$counts,NA), type='s', lwd=1.0, col=col)
        lines(h$breaks, c(NA,h$counts), type='h', lwd=1.0, col=col)
        lines(h$breaks, c(h$counts,NA), type='h',lwd=1.0, col=col)
        lines(h$breaks, rep(0,length(h$breaks)), type='S', col=col)
        invisible(h)
        
        dev.off()
        
        #######
        ##### Delays between mRNAs vs. proteins
        #######
        aa = nuclear.all[table.sx[,1], ]
        
        length(which(aa$qv.mrna<0.1))
        length(which(aa$qv<0.05 & aa$qv.mrna<=1))
        length(which(aa$qv<0.05 & aa$qv.mrna<0.1))
        length(which(aa$qv<0.05 & aa$q.total<=1))
        length(which(aa$qv<0.05 & aa$q.total<0.25))
        
        kk = which(aa$qv<0.05 & aa$qv.mrna<=1)
        nuclear.mrna = aa[kk,]
        nuclear.mrna$localization = table.sx$Localization[kk]
        nuclear.mrna$tfs = table.sx$TFs[kk]
        
        #### Compare phases of Nuclear and mRNA
        cutoff = 0.05
        cutoff.mrna = 0.05
        jj = which(nuclear.mrna$qv<cutoff & nuclear.mrna$qv.mrna<cutoff.mrna)
        print(length(jj))
        
        xx = data.frame(nuclear.mrna$phase.mrna[jj], nuclear.mrna$phase[jj], as.character(nuclear.mrna$localization[jj]), nuclear.mrna$Gene.names[jj], nuclear.mrna$tfs[jj], stringsAsFactors=FALSE)
        colnames(xx) = c('phase.mrna', 'phase.nuclear', 'localization', 'gene', 'tfs')
        
        kk = which(xx$localization=='Nucleus' | xx$localization=='Nucleus/Cytoplasm')
        pp = xx[kk,]
        o1 = order(pp[,1])
        pp = pp[o1,]
        circ.cor(pp[,1]/24*2*pi, pp[,2]/24*2*pi, test=TRUE)
        
        phyper(25, 67, (522-67),89,lower.tail=FALSE)
        genes.test = data.frame(unlist(strsplit(as.character(pp[,4]), ';')), stringsAsFactors=FALSE)
        genes.bg = data.frame(unlist(strsplit(as.character(aa$Majority.protein.IDs), ';')), stringsAsFactors=FALSE)
        #write.table(genes.test, file='Tables_DATA/R_R_genes.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        #write.table(genes.bg, file='Tables_DATA/R_R_genes_background_2.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        
        pdfname = paste('myplots/FIGURES/Phase_Correlation_mRNA_nuclear_both.pdf', sep='')
        pdf(pdfname, width=1.9, height=1.9)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        plot(pp[,1],c(nrow(pp):1), type='l', col='black',lwd=2.0, xlim=c(min(pp[,1]), (max(pp[,1])+24)), xlab=NA, ylab=NA, axes=FALSE)
        points(pp[,1]+24,c(nrow(pp):1), type='l', col='black', lwd=2.0)
        
        
        for(n in 1:nrow(pp))
        {
            if(pp[n, 3]=='Nucleus')
            {
                col1=rgb(0.2,1.0,0.2);
                pch = 19
                points(pp[n,2], c(nrow(pp):1)[n], col=col1, cex=0.5, pch=pch)
                points(pp[n,2]+24, c(nrow(pp):1)[n], col=col1,cex=0.5, pch=pch)
                
                
            }
            if(pp[n, 3]=='Nucleus/Cytoplasm')
            {
                #col=rgb(0.6,0,0.2)
                col2=rgb(0.1,0.6,0.1);
                pch = 19
                points(pp[n,2], c(nrow(pp):1)[n], col=col2, cex=0.5,pch=pch)
                points(pp[n,2]+24, c(nrow(pp):1)[n], col=col2,cex=0.5, pch=pch)
            }
            
        }
        axis(1,at=seq(0,48,by=12),cex.axis =1.0)
        axis(2,at=seq(0, 80, by=20), las=1,cex.axis = 1.0)
        box()
        #legend('topright', legend = c(' ',' ') , pch = c(1, 1), cex=0.7, col= c(col1, col2), border = NA, bty = 'n')
        
        dev.off()
        
        pp.diff = pp[,2] - pp[,1]
        kk = which(pp.diff<=(-12))
        pp.diff[kk] = pp.diff[kk]+24
        kk = which(pp.diff>12)
        pp.diff[kk] = pp.diff[kk]-24
        
        pp = cbind(pp, pp.diff)
        
        ii = which(pp[,3]=='Nucleus')
        jj = which(pp[,3]=='Nucleus/Cytoplasm')
        
        pdfname = paste('myplots/FIGURES/Phase_Delay_mRNA_nuclear.pdf', sep='')
        pdf(pdfname, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        hist(pp.diff[ii], breaks=seq(-12, 12, by=2), col=col1, xlab=NA, ylab=NA,main=NA, axes=FALSE)
        axis(1,at=seq(-12,12,by=6),cex.axis =1.0)
        axis(2,at=seq(0,16,by=4),las=1,cex.axis = 1.0)
        abline(v=median(pp.diff[ii]), col='darkgray', lwd=1.5)
        
        dev.off()
        
        pdfname = paste('myplots/FIGURES/Phase_Delay_mRNA_both.pdf', sep='')
        pdf(pdfname, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        hist(pp.diff[jj], breaks=seq(-12, 12, by=2), col=col2, xlab=NA, ylab=NA,main=NA, axes=FALSE)
        axis(1,at=seq(-12,12,by=6),cex.axis =1.0)
        axis(2,at=seq(0,8,by=4),las=1,cex.axis = 1.0)
        #box()
        abline(v=median(pp.diff[jj]), col='darkgray', lwd=1.5)
        dev.off()
        
        ks.test(pp.diff[ii], pp.diff[jj], alternative = c("two.sided"))
        
        ### delays of TFs
        pdfname = paste('myplots/FIGURES/Phase_Delay_mRNA_nuclear_TFs.pdf', sep='')
        pdf(pdfname, width=2.2, height=2.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        kk = which(pp$tfs==1)
        yy = pp[kk,]
        range(yy$pp.diff)
        o1 = order(yy$pp.diff)
        yy = yy[o1,]
        col1=rgb(0.2,1.0,0.2);
        col2=rgb(0.1,0.6,0.1);
        plot(yy$pp.diff, c(1:nrow(yy)), xlim=c(-15, 15), type='n', xlab=NA, ylab=NA,main=NA, axes=FALSE)
        abline(v=0, col='black', lwd=1.)
        abline(v=2, col='gray', lwd=1.5, lty=2)
        abline(v=-2, col='gray', lwd=1.5, lty=2)
        
        for(n in 1:nrow(yy))
        {
            if(yy$localization[n]=='Nucleus')
            {
                cols = col1
            }else{
                cols = col2
            }
            pos = 4
            if(yy$pp.diff[n]<0) pos=2
            points(c(yy$pp.diff[n], 0), c(n, n), type='l', col=cols, lwd=1.5)
            if(yy$gene[n]=='Trp53;Tp53')
            {
                text(yy$pp.diff[n], n, 'TP53', col='black', pos=pos, cex=0.6, offset=0.05)
            }else{
                text(yy$pp.diff[n], n, toupper(yy$gene[n]), col='black', pos=pos, cex=0.6, offset=0.05)
            }
            
        }
        axis(1,at=seq(-12,12,by=6),cex.axis = 1.0)
        
        #axis(2,at=seq(0,25,by=5),las=1,cex.axis = 0.8)
        
        dev.off()
        
        cutoff = 0.05
        cutoff.mrna = 0.05
        jj = which(nuclear.mrna$qv<cutoff & nuclear.mrna$qv.mrna<cutoff.mrna)
        print(length(jj))
        
        xx = data.frame(nuclear.mrna$phase.mrna[jj], nuclear.mrna$phase[jj], as.character(nuclear.mrna$localization[jj]), nuclear.mrna$Gene.names[jj], nuclear.mrna$tfs[jj], stringsAsFactors=FALSE)
        colnames(xx) = c('phase.mrna', 'phase.nuclear', 'localization', 'gene', 'tfs')
        
        
        pdfname = paste('myplots/FIGURES/Phase_Correlation_mRNA_Cyto_Cytoskeleton.pdf', sep='')
        pdf(pdfname, width=1.9, height=1.9)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        kk = which(xx$localization=='Cytoplasm' | xx$localization=='Cyto.Cytoskeleton')
        pp = xx[kk,]
        o1 = order(pp[,1])
        pp = pp[o1,]
        circ.cor(pp[,1]/24*2*pi, pp[,2]/24*2*pi, test=TRUE)
        
        plot(pp[,1],c(nrow(pp):1), type='l', col='black',lwd=2.0, xlim=c(min(pp[,1]), (max(pp[,1])+24)), xlab=NA, ylab=NA, axes=FALSE)
        points(pp[,1]+24,c(nrow(pp):1), type='l', col='black', lwd=2.0)
        
        for(n in 1:nrow(pp))
        {
            if(pp[n, 3]=='Cytoplasm')
            {
                col=rgb(1.0,0.4,0.5);
                pch = 19
                points(pp[n,2], c(nrow(pp):1)[n], col=col, cex=0.5,pch=pch)
                points(pp[n,2]+24, c(nrow(pp):1)[n], col=col,cex=0.5, pch=pch)
                
                
            }
            if(pp[n, 3]=='Cyto.Cytoskeleton')
            {
                col=rgb(0.6,0,0.2)
                pch = 19
                points(pp[n,2], c(nrow(pp):1)[n], col=col, cex=0.5,pch=pch)
                points(pp[n,2]+24, c(nrow(pp):1)[n], col=col,cex=0.5, pch=pch)
            }
            
        }
        axis(1,at=seq(0,48,by=12),cex.axis =1.0)
        axis(2,las=1,cex.axis = 1.0)
        box()
        
        dev.off()
        
        pp.diff = pp[,2] - pp[,1]
        kk = which(pp.diff<=(-12))
        pp.diff[kk] = pp.diff[kk]+24
        kk = which(pp.diff>12)
        pp.diff[kk] = pp.diff[kk]-24
        
        
        pdfname = paste('myplots/FIGURES/Phase_Delay_mRNA_Cyto_Cytoskeleton.pdf', sep='')
        pdf(pdfname, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        hist(pp.diff[ii], breaks=seq(-12, 12, by=2), col='gray', xlab=NA, ylab=NA,main=NA, axes=FALSE)
        axis(1,at=seq(-12,12,by=4),cex.axis =0.8)
        axis(2,at=seq(0,4,by=2),las=1,cex.axis = 0.8)
        
        dev.off()
        
        
        ##### Heatmap for proteins and mRNAs
        aa = nuclear.all[table.sx[,1], ]
        kk = which(aa$qv<0.05 & aa$qv.mrna<=1)
        nuclear.mrna = aa[kk,]
        nuclear.mrna$localization = table.sx$Localization[kk]
        nuclear.mrna$tfs = table.sx$TFs[kk]
        
        kk = which(aa$nb_timepoint.total>=8)
        nuclear.total = aa[kk,]
        nuclear.total$localization = table.sx$Localization[kk]
        nuclear.total$tfs = table.sx$TFs[kk]
        
        require('CircStats')
        locas = c('mix', 'cyto')
        
        for(loca in locas)
        {
            ### compare mRNA and nuclear proteins
            loca = 'mix'
            if(loca =='mix')
            {
                jj = which(nuclear.mrna$localization=='Nucleus'|nuclear.mrna$localization=='Nucleus/Cytoplasm')
            }else{
                jj = which(nuclear.mrna$localization=='Cytoplasm'|nuclear.mrna$localization=='Cyto.Cytoskeleton')
            }
            
            xx = nuclear.mrna[jj,]
            #### Heatmaps for mRNAs and rhythmic nuclear proteins
            library("gplots")
            cutoff = 0.05 ## FDR for rhythmic nuclear proteins
            cutoff.mrna = 0.05 ## FDR for rhythmic mRNAs
            
            dim(xx)
            kk = which(xx$qv.mrna<cutoff.mrna)
            ii = which(xx$qv.mrna>=cutoff.mrna)
            print(length(kk))
            print(length(ii))
            
            o1 = order(xx$phase[kk])
            kk = kk[o1]
            o2 = order(xx$phase[ii])
            ii = ii[o2]
            
            
            x = xx[kk, 1:16]
            x = rbind(x,rep(NA,ncol(x)))
            x = rbind(x,rep(NA,ncol(x)))
            x = rbind(x,xx[ii, 1:16])
            x = t(apply(x, 1, standadize.nona))
            #x = rbind(x,nuclear.mrna[ii, grep('mRNA.ZT', colnames(nuclear.mrna))])
            
            y=xx[kk,grep('mRNA.ZT', colnames(xx))]
            y = rbind(y,rep(NA,ncol(y)))
            y = rbind(y,rep(NA,ncol(y)))
            y = rbind(y, xx[ii, grep('mRNA.ZT', colnames(xx))])
            y = t(apply(y, 1, standadize.nona))
            
            x = as.matrix(x)
            y = as.matrix(y)
            #x = x[, c(1:7)]
            
            pdfname = paste('myplots/FIGURES/Heatmap_mRNA_', loca, '_GO.pdf', sep='')
            if(loca=='mix')
            {
                pdf(pdfname, width=0.7, height=2)
            }
            if(loca=='cyto')
            {
                pdf(pdfname, width=0.7, height=0.7)
            }
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0,ColSideColors=c(rep('gray70',6),rep('black',6),rep('gray70',6),rep('black',6),'white'),margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,2.5),
            #		  labRow=NA,cexRow=0.5,xlab=NA, key=TRUE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
            
            heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA,  na.rm = TRUE, labCol=NA, keysize=1.0, margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05, 3.0),
            labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
            
            dev.off()
            
            pdfname = paste('myplots/FIGURES/Heatmap_Nucelar_', loca, '_GO.pdf', sep='')
            if(loca=='mix')
            {
                pdf(pdfname, width=0.7, height=2)
            }
            if(loca=='cyto')
            {
                pdf(pdfname, width=0.7, height=0.7)
            }
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0,ColSideColors=c(rep('gray70',4),rep('black',4),rep('gray70',4),rep('black',4)),margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,2.5),
            #		  labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
            heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0,margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
            labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
            
            dev.off()
            
            #### Compare total and nuclear proteins
            loca = 'cyto'
            if(loca =='mix')
            {
                jj = which(nuclear.total$localization=='Nucleus'|nuclear.total$localization=='Nucleus/Cytoplasm')
            }else{
                jj = which(nuclear.total$localization=='Cytoplasm'|nuclear.total$localization=='Cyto.Cytoskeleton')
            }
            
            xx = nuclear.total[jj,]
            #### Heatmaps for mRNAs and rhythmic nuclear proteins
            library("gplots")
            cutoff = 0.05 ## FDR for rhythmic nuclear proteins
            cutoff.total = 0.25 ## FDR for rhythmic mRNAs
            
            kk = which(xx$q.total<cutoff.total)
            ii = which(xx$q.total>=cutoff.total)
            print(nrow(xx))
            print(length(kk))
            print(length(ii))
            
            ### let's see if we need compare the phases for total and nuclear proteins or not
            ssp = data.frame(xx$Gene.names[kk], xx$phase[kk], xx$phase.total[kk], xx$localization[kk], stringsAsFactors=FALSE)
            
            phases = ssp[,2]
            pdf('myplots/FIGURES/Phase_distribution_nuclear_vs_total_cyto_cytoskeleton_1.pdf', width=1.6, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            circular_phase24H_histogram(phases, col=rgb(0.6,0.,0.2), cex.axis=0.5, cex.lab=0.01, lwd=0.5)
            
            dev.off()
            
            phases = ssp[,3]
            pdf('myplots/FIGURES/Phase_distribution_nuclear_vs_total_cyto_cytoskeleton_2.pdf', width=1.6, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            circular_phase24H_histogram(phases, col=rgb(0.2,0.2,0.5), cex.axis=0.5, cex.lab=0.01, lwd=0.5)
            
            dev.off()

            
            o1 = order(xx$phase[kk])
            kk = kk[o1]
            o2 = order(xx$phase[ii])
            ii = ii[o2]
            
            
            x = xx[kk, 1:16]
            x = rbind(x,rep(NA,ncol(x)))
            x = rbind(x,rep(NA,ncol(x)))
            x = rbind(x,xx[ii, 1:16])
            x = t(apply(x, 1, standadize.nona))
            #x = rbind(x,nuclear.mrna[ii, grep('mRNA.ZT', colnames(nuclear.mrna))])
            
            y=xx[kk,grep('.total', colnames(xx))[1:16]]
            y = rbind(y,rep(NA,ncol(y)))
            y = rbind(y,rep(NA,ncol(y)))
            y = rbind(y, xx[ii, grep('.total', colnames(xx))[1:16]])
            y = t(apply(y, 1, standadize.nona))
            
            x = as.matrix(x)
            y = as.matrix(y)
            #x = x[, c(1:7)]
            
            pdfname = paste('myplots/FIGURES/Heatmap_Total_prot_', loca, '_GO.pdf', sep='')
            pdf(pdfname, width=0.5, height=1.7)
            
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0,ColSideColors=c(rep('gray70',6),rep('black',6),rep('gray70',6),rep('black',6),'white'),margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,2.5),
            #		  labRow=NA,cexRow=0.5,xlab=NA, key=TRUE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
            
            heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA,  na.rm = TRUE, labCol=NA, keysize=1.0, margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05, 3.0),
            labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
            
            dev.off()
            
            pdfname = paste('myplots/FIGURES/Heatmap_Total_Nucelar_prot_', loca, '_GO.pdf', sep='')
            pdf(pdfname, width=0.5, height=1.7)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0,ColSideColors=c(rep('gray70',4),rep('black',4),rep('gray70',4),rep('black',4)),margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,2.5),
            #		  labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
            heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0,margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
            labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
            
            dev.off()
            
            
        }
        
    }
    
    #######
    ##### functional analysis of Nuclear, Nuclear/Cytoplasmic, Cytoplasmic proteins
    #######
    functional.analysis = FALSE
    if(functional.analysis)
    {
        jj = which(table.sx$Localization=='Nucleus')
        genes = table.sx$Gene.names[jj]
        genes = unlist(strsplit(as.character(genes), ';'))
        genes = unlist(strsplit(as.character(genes), '/'))
        write.table(genes, file='Tables_DATA/Table_SX_Nucleus.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        
        jj = which(table.sx$Localization=='Nucleus/Cytoplasm')
        genes = table.sx$Gene.names[jj]
        genes = unlist(strsplit(as.character(genes), ';'))
        genes = unlist(strsplit(as.character(genes), '/'))
        write.table(genes, file='Tables_DATA/Table_SX_Nucleus_Cytoplam.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        
        jj = which(table.sx$Localization=='Cytoplasm' |table.sx$Localization=='cytoplasm')
        genes = table.sx$Gene.names[jj]
        genes = unlist(strsplit(as.character(genes), ';'))
        genes = unlist(strsplit(as.character(genes), '/'))
        write.table(genes, file='Tables_DATA/Table_SX_Cytoplam.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
        
        locas = c('Nucleus', 'Nucleus_Cytoplasm', 'Cytoplasm')
        sel = c(1:12)
        for(loca in locas)
        {
            file = paste('Tables_DATA/David_Analysis_', loca, '.txt', sep='')
            aa = read.delim(file, sep='\t', header=FALSE)
            kk = grep('Annotation Cluster', aa[,1])
            enrich = as.numeric(unlist(strsplit(as.character(aa[kk,2]), ':'))[c(1:length(kk))*2])
            jj = which(enrich>2)
            #print(length(jj))
            enrich = enrich[sel]
            ffn = c()
            for(n in 1:length(sel))
            {
                bb = aa[c(kk[n]:kk[n+1]),]
                jj = grep('SP_PIR_KEYWORDS', bb[,1])
                if(length(jj)>0)
                {
                    jj = jj[1];
                    ffn = c(ffn, as.character(bb[jj,2]))
                }else
                {
                    jj = grep('GOTERM', bb[,1])
                    if(length(jj)>0)
                    {
                        if(n==7 & loca=='Nucleus') jj=jj[2]
                        jj = jj[1];
                        ffn = c(ffn, as.character(unlist(strsplit(as.character(bb[jj, 2]), '~'))[2]))
                    }else
                    {
                        ff = as.character(unlist(strsplit(as.character(bb[3, 2]), ':'))[2])
                        if(length(ff)==0)
                        {
                            print(n);
                            print(loca)
                            print('not found')
                            ffn = c(ffn, NA)
                        }else{ 
                            ffn = c(ffn, ff)
                        }
                    }
                }
            }
            
            o1 = order(enrich)
            enrich = enrich[o1]
            ffn = ffn[o1]
            
            pdfname = paste('myplots/FIGURES/Function_analysis_David_', loca, '.pdf', sep='')
            pdf(pdfname, width=2.2, height=1.5)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,7.,0.5,0.5)+0.1, tcl = -0.3)
            
            barplot(enrich, main=NA, col=c("gray"),names.arg=ffn, beside=TRUE, horiz = TRUE, space=0.25, las=1, cex.lab=0.1,cex.names=0.5)
            
            dev.off()
            
        }
    }
    
    #############
    #### Examples of proteins with different localizations
    #############
    hmr = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/DATA_WB.csv', sep=',', header=TRUE)
    colnames(hmr) = toupper(colnames(hmr))
    
    require('plotrix')
    
    examples1 = c('Gsk3a', 'Nr3c1', 'Anxa2','Lcp1',  'Acaca','Fasn', 'Alb')
    cols = 'plum'
    
    for(gene in examples1)
    {
        pdf.name = paste("myplots/FIGURES/Fig_2_Examples_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.9, height=1.6)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        n = which(nuclear$Gene.names==gene)
        y0 = as.numeric(nuclear[n, c(1:16)])
        y0 = 2^y0
        y0 = y0/mean(y0, na.rm=TRUE)
        
        kk = grep(toupper(gene), colnames(hmr))
        
        y1 = rbind(as.numeric(unlist(hmr[c(1:16), kk[1]])), as.numeric(unlist(hmr[c(1:16), kk[2]])))
        #y1 = rbind(y1[c(1:8)], y1[c(9:16)], y1[17:24], y1[25:32])
        for(ii in 1:2)
        {
            y1[ii, ] = y1[ii, ]/mean(y1[ii, ], na.rm=TRUE)
        }
        
        averg1 = c()
        err1 = c()
        for(ii in 1:16)
        {
            xx = y1[,ii]
            #xx = xx/mean(xx[which(!is.na(xx)==TRUE)])
            ll = length(which(is.na(xx)==TRUE))
            if(ll==0){averg1 = c(averg1, mean(xx)); err1=c(err1, sd(xx)/sqrt(2))}
            if(ll==1){averg1 = c(averg1, mean(xx[which(!is.na(xx)==TRUE)])); err1=c(err1, sd(xx[which(!is.na(xx)==TRUE)])/sqrt(length(xx[which(!is.na(xx)==TRUE)])))}
            if(ll==2){averg1 = c(averg1, NA); err1=c(err1, NA)}
            
        }
        
        time = c(0:15)*3
        lims = range(c(y0, averg1-err1, averg1+err1), na.rm=TRUE)
        
        #lims = range(c(y0[which(!is.na(y0)==TRUE)]))
        pch = 16
        plot(time, y0, type='l', ylim=lims, xlim=c(0, 45), col='darkblue', lwd=1., pch=pch, main=toupper(gene), cex.main=0.8, axes=FALSE, xlab=NA, ylab=NA)
        points(time, y0, type='p', cex=0.7, col='darkblue', pch=16)
        points(time, averg1, type='l', lwd=1, col=cols, pch=15)
        points(time, averg1, type='p', cex=0.8, col=cols, pch=15)
        #arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col='darkblue', lwd=1.5)
        arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3, col=cols, lwd=1)
        
        #axis(1,at=12*c(0:4),cex.axis =1.2)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
        by = 0.1
        print(gene)
        print(lims)
        
        cex.axis = 0.8
        #axis(2,at = seq(lims[1], lims[2],by=0.3),las=1,cex.axis = 1.2)
        if(gene=='Gsk3a') axis(2,at = c(0.6, 1, 1.4),las=1,cex.axis = cex.axis)
        if(gene=='Nr3c1') axis(2,at = c(0.6, 1, 1.4),las=1,cex.axis = cex.axis)
        if(gene=='Anxa2') axis(2,at = c(0.5, 1, 1.5),las=1,cex.axis = cex.axis)
        if(gene=='Lcp1') axis(2,at = c(0.5, 1, 1.5, 2.0),las=1,cex.axis = cex.axis)
        if(gene=='Acaca') axis(2,at = c(0.5, 1, 1.5, 2.0, 2.5, 3.0),las=1,cex.axis = cex.axis)
        if(gene=='Fasn') axis(2,at = c(0.5, 1, 1.5, 2.0, 2.5, 3.0),las=1,cex.axis = cex.axis)
        if(gene=='Alb') axis(2,at = c(0.5, 1, 1.5, 2.0, 2.5, 3.0),las=1,cex.axis = cex.axis)
        
        box()
        
        abline(h=1,lty=2,lwd=1.5, col="darkgray")
        #if(gene=='Gsk3a' | gene=='Alb')
        #{
        #    legend(36, (lims[2]), c('MS','WB'),pch=c(16,15),lty=c(1,1),col=c("darkblue",cols), cex=0.8, pt.lwd=1.5, pt.cex=1.2,bty='n')
        #}
        dev.off()
        
    }
    
    ######
    ### confocal signals for Fasn and Acaca
    ######
    conforcal = read.csv('Tables_DATA/DATA_CONFOCAL_FASN.csv')
    ss = c()
    for(n in c(0:15)*3)
    {
        kk = which(conforcal[,5]==n)
        test = conforcal[kk, ]
        index.n = unique(test[,3])
        nn = c()
        for(index in index.n)
        {
            jj = which(test[,3]==index)
            nn = c(nn, mean(log2(as.numeric(test[jj, 4]))), na.rm=TRUE)
        }
        ss = c(ss, mean(nn, na.rm=TRUE))
    }
    #ss = 2^ss
    ### centralizing two samples
    #ss[1:8] = ss[1:8]-mean(ss[1:8])
    #ss[9:16] = ss[9:16]-mean(ss[9:16])+mean(ss[1:8])
    
    pdf('myplots/FIGURES/FASN_conforcal_signals.pdf', width=1.8, height=1.6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    source('functions_nuclear.R')
    #averg0 = mean.err(keep$cell.size.all.mean, period=48)[1,]
    #err0 = mean.err(keep$cell.size.all.mean, period=48)[2,]
    averg = mean.err(ss, period=24)[1,]
    err = mean.err(ss, period=24)[2,]
    
    lims = range(c(averg+err, averg-err));
    pch=19
    lwd = 1.
    time = c(0:7)*3
    xlim = c(0, 24)
    col = 'black'
    cex = 0.6
    
    plot(time, averg, type='p', ylim=lims, xlim=xlim, col=col, lwd=1., cex=cex, lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1., cex=cex)
    col='black'
    
    tt = seq(0, 21, by=0.2)
    res = f24_R2_alt2(as.numeric(ss), t=c(0:15)*3)
    print(res)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.5)
    
    #cex = 0.8;
    axis(1,at=6*c(0:4),cex.axis =1.0)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    axis(2, at= c(-0.6, -0.4, -0.2, 0, 0.2), las=1,cex.axis = 1.0)
    box()
    #legend(7, 4500, c('mononu', 'binu'), pch=c(1, 0), lty=c(1,2),col=c('blue', 'red'), lwd=1.5, cex=0.9, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()

    
    conforcal = read.csv('Tables_DATA/DATA_CONFOCAL_ACACA.csv')
    ss = c()
    for(n in c(0:15)*3)
    {
        kk = which(conforcal[,5]==n)
        test = conforcal[kk, ]
        index.n = unique(test[,3])
        nn = c()
        for(index in index.n)
        {
            jj = which(test[,3]==index)
            nn = c(nn, mean(log2(as.numeric(test[jj, 4]))), na.rm=TRUE)
        }
        ss = c(ss, mean(nn, na.rm=TRUE))
    }
    #ss = ss^2
    pdf('myplots/FIGURES/ACACA_conforcal_signals.pdf', width=1.8, height=1.6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    source('functions_nuclear.R')
    #averg0 = mean.err(keep$cell.size.all.mean, period=48)[1,]
    #err0 = mean.err(keep$cell.size.all.mean, period=48)[2,]
    averg = mean.err(ss, period=24)[1,]
    err = mean.err(ss, period=24)[2,]
    
    lims = range(c(averg+err, averg-err));
    pch=19
    lwd = 1.
    time = c(0:7)*3
    xlim = c(0, 24)
    col = 'black'
    cex = 0.6
    plot(time, averg, type='p', ylim=lims, xlim=xlim, col=col, lwd=1., cex=cex, lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1., cex=cex)
    col='black'

    tt = seq(0, 24, by=0.2)
    res = f24_R2_alt2(as.numeric(ss), t=c(0:15)*3)
    print(res)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=2.0)
    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    axis(2, at= c(0.2, 0.4, 0.6), las=1,cex.axis = 1.0)
    box()
    #legend(7, 4500, c('mononu', 'binu'), pch=c(1, 0), lty=c(1,2),col=c('blue', 'red'), lwd=1.5, cex=0.9, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()
    
}

#################
####### FIGURE 3 (Main + SUpplementary)
#################
FIGURE_3 = TRUE
if(FIGURE_3)
{
    nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
    nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')
    #load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Table_Nuclear_mRNA_Total_Nascent.Rdata')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    load(file='Rdata/Anotation_protein_complex_Corum_mannual_detected_nuleus.Rdata')
    load(file='Tables_DATA/Table_Nuclear_total_mRNA_pol2.Rdata')
    load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/Annotation_Used_All_V2.Rdata')
    load(file='Tables_DATA/Table_mRNA_CRY_WT_NRF.Rdata')

    annot = annot.corum
    
    ii = which(!is.na(annot$localization) & annot$localization!='nucleus')
    length(ii)
    kk = which(annot$nb.subunits>0)
    annot = annot[kk,]
    
    ############
    #### Our detected nuclear proteins are mainly detected in the nuclear protein complexes
    ###########
    nucleus = read.delim('Tables_DATA/Annotation_pure_nucleus_proteins_v3.txt', sep='\t', header=FALSE,as.is=TRUE)
    mix = read.delim('Tables_DATA/Annotation_mix_cytoplasmic_nucleus_proteins_v3.txt', sep='\t', header=FALSE, as.is=TRUE)
    #cyto = read.delim('Tables_DATA/Annotation_pure_cytoplasmic_proteins_v2.txt', sep='\t', header=TRUE)
    nc = c(nucleus[,1], mix[,1])
    nc = unique(nc)
    
    mapping.all = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Manually_curation_TFs_pools.txt', header=TRUE, sep='\t')
    poly = read.table('/Users/jiwang/RNA_seq_Data/Total-RNA-seq/mRNAs_pre_mRNAs_PolyA_RNA_seq_Gachon.txt', header=TRUE, sep='\t') ### this polyA RNA-seq data needs to be checked
    ensgene = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/ensgenes_mapping.txt')
    cutoff.expressed = 0
    gene.expressed = unique(c(as.character(poly$gene[which(poly$mean.ex>=cutoff.expressed)]), as.character(nuclear.names[,3])))
    ### expressed genes annotated as nuclear or nuclear/cytoplasmic
    nc.expressed = intersect(nc, gene.expressed)
    
    nc.complex = annot$subunits
    nc.complex = unlist(strsplit(nc.complex, ','))
    nc.complex = unique(nc.complex)
    nc.expressed.complex = intersect(nc.complex, nc.expressed)
    nc.expressed.solo = setdiff(nc.expressed, nc.complex)
    
    length(which(!is.na(match(nc.expressed.solo, nuclear.names[,3]))))/length(nc.expressed.solo)
    length(which(!is.na(match(nc.expressed.complex, nuclear.names[,3]))))/length(nc.expressed.complex)
    
    nb.subunits = unique(annot$nb.subunits)
    nb.subunits = nb.subunits[order(nb.subunits)]
    
    keep = c()
    
    for(n in 1:nrow(annot))
    #for(nb in nb.subunits)
    {
        xx = annot$subunits[n]
        xx = unlist(strsplit(xx, ','))
        xx = unique(xx)
        keep = c(keep, length(intersect(xx, nc.expressed))/length(xx))
        
    }
    #keep = data.frame(keep0, keep1, keep2, keep, stringsAsFactors=FALSE)
    
    pdf.name = paste("myplots/FIGURES/Enrichment_nuclear_prots_PCs.pdf")
    pdf(pdf.name, width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    kk = which(keep>0.25)
    boxplot(annot$percent.detected[kk] ~ annot$nb.subunits[kk], ylim=c(0, 1.0), axes=FALSE, pars = list(boxwex = 0.7, staplewex = 0.4, outwex = 0.4), cex.axis=1.)
    
    index = unique(annot$nb.subunits[kk])
    index = index[order(index)]
    labels = c(1, 10, 20, 30, 80)
    xx = match(labels, index)
    axis(1,at=xx,labels=labels, cex.axis =1.)
    axis(2,at=seq(0, 1, by=0.2), las=1,cex.axis = 1.)
    #box()
    abline(h=length(which(!is.na(match(nc.expressed.solo, nuclear.names[,3]))))/length(nc.expressed.solo), lwd=1.5, col='gray')
    
    box()
    
    dev.off()
    
    #########
    #### SVD analysis to determine expressed and rhythmic complexes
    #########
    #### Examples of PER and 20S
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
    load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v2.Rdata' )
    aa = res
    kk = which(aa$pval.svd<0.05)
    aa = aa[kk,]
    #Plots.complexes.details(aa, c(1:nrow(aa)), 'myplots/All_PCs_Corum_subunits_Selected_SVD_pval_0.05_details.pdf')
    
    indexx = c(20, 4)
    pc.names = c('26S-Proteasome','PER complex')
    pc.cols = c(rep('black', 2))
    
    for(n in 1:length(indexx))
    {
        kk = indexx[n]
        name = pc.names[n]
        index = aa$index.detected[kk]
        index = as.numeric(unlist(strsplit(as.character(index), ',')))
        genes = unlist(strsplit(as.character(aa$subunits.detected[kk]), ','))
        index = index[which(nuclear$nb.timepoints[index]>=8)]

        pdf.name = paste("myplots/FIGURES/SVD_examples_sup_profiles_", name, ".pdf", sep='')
        pdf(pdf.name, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        if(length(index)>1)
        {
            test = as.matrix(nuclear[index,c(1:16)])
            test =apply(test, 1, standadization.nona.impute)
            #test =apply(test, 1, standadize.nona)
            test = t(test)
            pvals = signif(nuclear$pval[index],d=2)
            phases = signif(nuclear$phase[index],d=2)
            amp = signif(nuclear$relamp[index],d=1)
            o1 = order(pvals)
            pvals = pvals[o1]
            test = test[o1,]
            phases = phases[o1]
            index = index[o1]
            genes = genes[o1]
            ylim = range(test, na.rm=TRUE)
            xlim = c(0,48)
            
            plot(1,1,type = 'n', xlab = NA, ylab = NA, xlim = xlim, ylim = ylim, main=paste(pc.names[n]), axes=FALSE, cex.main=0.8)
            abline(h=0, col='gray', lty=2, lwd=1.5)
            
            #text.col = c(rep('cornflowerblue', 5),  rep('darkorchid',5))
            #text.col = rep('blue', length(genes))
            text.col = 'black'
            for(i in 1:nrow(test))
            {
                i.rel = i;
                type.plot='b';
                col= 'gray'
                
                points(c(0:15)*3, test[i,], type='l', lwd=1., col = col , lty = 1)
                #legend(x = 41,y = (ylim[2]-(i-1)*0.2), legend = paste(genes[i]), text.col=text.col, bty = 'n', cex=0.75)
            }
            median = apply(test, 2, mean)
            points(c(0:15)*3, median, type='l', col=pc.cols[n], lwd=2.0)
            
            axis(1,at=seq(0, 48, by=12),cex.axis =1.)
            if(n==1) axis(2,at = seq(-2, 2, by=1),las=1,cex.axis = 1.)
            if(n==2) axis(2,at = seq(-2, 2, by=1),las=1,cex.axis = 1.)
            box()
        }
        dev.off()
        
        #### SVD analysis
        bg.cutoff = aa$bg.cutoff[kk]
        test = matrix(data = NA, nrow = length(index), ncol = 16)
        for(mm in 1:length(index))
        test[mm, ] = standadization.nona.impute(as.numeric(nuclear[index[mm], c(1:16)]))
        ss = svd(test)
        L = length(ss$d)
        pl = ss$d^2/sum(ss$d^2)
        
        u1 = ss$u[,1]
        u2 = ss$u[,2]
        v1 = ss$v[,1]
        v2 = ss$v[,2]
        nb.positive = length(which(u1>=0))
        nb.negative = length(which(u1<0))
        if(nb.positive<nb.negative){u1 = -u1; v1 = -v1;}
        
        res1 = pl[1]*length(which(u1>=0))/length(u1);
        stat1 = f24_R2_alt2(v1, t=c(0:15)*3)
        
        pdf.name = paste("myplots/FIGURES/SVD_examples_sup_distribution_components_", name, ".pdf", sep='')
        pdf(pdf.name, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        barplot(pl, col='darkgray', ylim=c(0,1), axes=FALSE)
        abline(h=bg.cutoff, col='darkgray',lwd=1.5)
        
        #if(name=='26S-Proteasome') axis(1,at=c(1, 5, 10, 14),cex.axis =1.)
        #if(name=='PER complex') axis(1,at=seq(1, 7, by=2),cex.axis =1.)
        
        axis(2,at = seq(0, 1, by=0.2),las=1,cex.axis = 1.)
        box()
        
        dev.off()
        
        pdf.name = paste("myplots/FIGURES/SVD_examples_sup_1st_components_", name, ".pdf", sep='')
        pdf(pdf.name, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        plot(c(0:15)*3, apply(ss$u, 2, mean)[1]*ss$v[,1], type='l',xlim=c(0, 48), lwd=1., cex=0.7, col='black',xlab=NA,ylab=NA,main=NA, axes=FALSE)
        points(c(0:15)*3, apply(ss$u, 2, mean)[1]*ss$v[,1], type='p', cex=0.7, col='black')
        abline(h=0, lty=2, col='darkgray',lwd=1.5)
        axis(1,at=seq(0, 48, by=12),cex.axis =1.)
        if(name=='26S-Proteasome') axis(2,at = c(-0.1, 0),las=1,cex.axis = 1.)
        if(name=='PER complex') axis(2,at = c(-0.1, 0, 0.1),las=1,cex.axis = 1.)
        box()
        dev.off()
						
    }
    
    #### SVD analysis for all protein complexes and identified rhythmic protein complexes
    filename = 'Rdata/Analysis_res_Protein_Complexes_SVD_Rhythmic_FDR_all_v3.Rdata'
    load(file=filename)
    pdf.name = paste("myplots/FIGURES/SVD_Complexes_v2.pdf")
    pdf(pdf.name, width=2.5, height=2.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    kk = which(res$qv.p1<0.05 & res$pval.svd<0.05)
    jj = which(res$qv.p1>=0.05 & res$pval.svd<0.05)
    ii = which(res$pval.svd>=0.05)
    nb.subunits = res$nb.subunits.quantified
    nb.subunits = unique(nb.subunits)
    nb.subunits = nb.subunits[order(nb.subunits)]
    bg.cutoff = res$bg.cutoff[match(nb.subunits, res$nb.subunits.quantified)]
    cex = 0.7
    plot(res$nb.subunits.quantified[kk], res$p1[kk], col='darkred',log='x', cex=cex, pch=1, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    points(res$nb.subunits.quantified[jj], res$p1[jj], col='blue', pch=1,cex=cex)
    points(res$nb.subunits.quantified[ii], res$p1[ii], col='black', pch=1,cex=0.4)
    points(nb.subunits, bg.cutoff, type='l', col='gray', lwd=2.0)
				
    axis(1,at=c(2, 5, 10, 20, 50, 100),cex.axis =1)
    axis(2,at = seq(0,1, by=0.2),las=1,cex.axis = 1)
    box()
    
    dev.off()
    
    Table.Sup = FALSE
    if(Table.Sup)
    {
        load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v4.Rdata')
        
        kk = which(res$pval.svd<0.05)
        xx = res[kk, ]
        yy = xx[, c(1,2,6,7,8,10,11, 14, 17, 18, 27, 26)]
        
        write.table(yy, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S3_nuclear_proteins_complexes.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        xx = res.c[, -c(5, 6)]
       
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S3_nuclear_proteins_complexes_rhythmic.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
    }

    
    #############
    ###### Phase Distribution of Rhythmic protein complexes
    #############
    load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v4.Rdata')
    res.c = res.c[, -c(5,6)]
    
    jj = which(res.c$amp.pc>1.5)
    
    kk = which(res.c$amp.pc<=2.1)
    res.c = res.c[kk,]
    #write.table(res.c, file='Tables_DATA/Phase_Amps_Function_Rhythmic_PCs.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    o2 = order(res.c$phase.pc)
    phase.pc = res.c$phase.pc[o2]
    names.pc = res.c$names.pc[o2]
    amp.pc = res.c$amp.pc[o2]
    #xx = res.c[o2,]
    #xx = xx[o2, ]
    
    Adjust.phase = FALSE
    if(Adjust.phase)
    {
        xx = phase.pc
        yy = phase.pc
        cutoff.adjust = 0.2
        for(n in 1:(length(xx)-1))
        {
            diff = (xx[n+1]-xx[n]);
            if(diff>=24) diff = diff-24;
            if(diff<cutoff.adjust) xx[n+1] = (xx[n]+cutoff.adjust)%%24;
            #if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
        }
        phase.pc = xx
        
    }else{
        length(which(phase.pc>=4.0 & phase.pc<=12))
        length(which(phase.pc<=2.5 | phase.pc>=19))
        pp.pc = phase.pc
        cutoff.adjust = 0.3
        kk = which(phase.pc>=4.0 & phase.pc<=12)
        for(n in 1:(length(kk)-1))
        {
            cat(n, '\n')
            diff = (phase.pc[kk[n+1]]-phase.pc[kk[n]]);
            if(diff>=24) diff = diff-24;
            if(diff<cutoff.adjust) pp.pc[kk[n+1]] = (pp.pc[kk[n]]+cutoff.adjust)%%24;
            #pp.pc[kk] = 4.0 + (12-4.0)/(length(kk)-1)*c(0:c(length(kk)-1))
        }
        
        kk = which(phase.pc<=2.5 | phase.pc>=19)
        #kk = c(kk[-c(1:15)], kk[c(1:15)])
        pp.pc[kk] = 19 + (24+2.5-19)/(length(kk)-1)*c(0:c(length(kk)-1))
        pp.pc = pp.pc%%24
        phase.pc = pp.pc
    }
    xx = res.c[o2,]
    tf.amp = 10
    #Plots.complexes(kk, 'myplots/PCs_50_percentage_rhythmic_subunits.pdf')
    tf.a = tf.amp*cos(2*pi/24*phase.pc)
    tf.b = tf.amp*sin(2*pi/24*phase.pc)
    CC= (tf.a -1i*tf.b) * exp(1i*pi/2)
    tf.a = Re(CC)
    tf.b = Im(CC)
    
    tf.amp = tf.amp + 0.4
    tf.aa = tf.amp*cos(2*pi/24*phase.pc)
    tf.bb = tf.amp*sin(2*pi/24*phase.pc)
    CC= (tf.aa -1i*tf.bb) * exp(1i*pi/2)
    tf.aa = Re(CC)
    tf.bb = Im(CC)
    
    rmm=max(abs(c(tf.aa,tf.bb)))*2
    rr=c(-rmm,rmm)
    
    bg.cols = c('red','darkturquoise', 'darkgreen', 'orange', 'lightslateblue', 'yellow', 'magenta', 'black')
    
    pdf('myplots/FIGURES/Phases_rhythmic_PCs_Corum_manual_v2.pdf', width = 6.2, height = 6.2)
    #par(cex = 0.5, las = 1, mgp = c(0,0,0), mar = c(0.1,0.1,0.1,0.1)+0.1, tcl = -0.3)
    #par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    plot(tf.a, tf.b, main=NA, type='n', xlim=rr, ylim=rr, axes=F, xlab='', ylab='', lwd=1.0, pch=1,cex=0.7)
    #dist = 0.2
    #arrows((tf.cycle)*(-1), 0, (tf.cycle), 0,  col='gray30', lty=4, lwd=1.5, length=0.0);
    #arrows(0, (tf.cycle)*(-1), 0, (tf.cycle), col='gray30', lty=1, lwd=1.5, length=0.0);
    
    #abline(v=0,h=0, col='darkgray',lwd=1.5)
    phi=seq(0,2*pi,len=1000)
    dist = 0.75
    lines((tf.amp-dist)*cos(phi), (tf.amp-dist)*sin(phi), col='darkgray', lwd=2.0)
    tf.cycle = tf.amp - dist
    
    lines(tf.cycle*4/5*cos(phi), 4/5*tf.cycle*sin(phi), col='darkgray', lwd=1.0,lty=2)
    lines(tf.cycle*3/5*cos(phi), 3/5*tf.cycle*sin(phi), col='darkgray', lwd=1.,lty=2)
    lines(tf.cycle*2/5*cos(phi), 2/5*tf.cycle*sin(phi), col='darkgray', lwd=1.,lty=2)
    lines(tf.cycle*1/5*cos(phi), 1/5*tf.cycle*sin(phi), col='darkgray', lwd=1.,lty=2)
    
    arrows((tf.amp-dist)*(-1), 0, (tf.amp-dist), 0,  col='darkgray', lty=2, lwd=1.0, length=0.0);
    arrows(0, (tf.amp-dist)*(-1), 0, (tf.amp-dist), col='darkgray', lty=2, lwd=1.0, length=0.0);
    
    for(n in 1:length(tf.a))
    {
        ii = n
        
        if(phase.pc[n]<=6) srt = 90-phase.pc[n]/24*360;
        if(phase.pc[n]>6 & phase.pc[n]<=12) srt = 90-phase.pc[n]/24*360;
        if(phase.pc[n]>12 & phase.pc[n]<=18) srt = 270-phase.pc[n]/24*360;
        if(phase.pc[n]>=18) srt = 270-phase.pc[n]/24*360;
        
        col='black'
        
        #if(!is.na(xx$tfs.pc[ii])==TRUE & is.na(xx$cors.pc[ii])==TRUE) col='magenta';
        #if(is.na(xx$tfs.pc[ii])==TRUE & !is.na(xx$cors.pc[ii])==TRUE) col='darkblue';
        #if(!is.na(xx$tfs.pc[ii])==TRUE & !is.na(xx$cors.pc[ii])==TRUE) col='magenta';
        
        if(phase.pc[n]<=12.0) {
            adj = 0;
        }else{
            adj = 1;
        }
        bg = 'black'
        if(xx$functions[n]=='Transcriptional regulation') bg = bg.cols[1]
        if(xx$functions[n]=='Cytoskeleton') bg = bg.cols[2]
        if(xx$functions[n]=='Chaperone') bg = bg.cols[3]
        if(xx$functions[n]=='Cell cycle') bg = bg.cols[4]
        if(xx$functions[n]=='DNA repair') bg = bg.cols[5]
        if(xx$functions[n]=='Proteolisis') bg = bg.cols[6]
        if(xx$functions[n]=='Transport') bg = bg.cols[7]
        
        points(tf.a[n], tf.b[n], type='p', pch=21, col=col, bg=bg, cex=1.2, lwd=1.0)
        text(tf.aa[n], tf.bb[n], names.pc[n], cex=0.7, col=col, srt=srt, adj=adj)
        
        amp = (1-(1-log2(amp.pc[n]))/5)*tf.cycle
        amp.a = amp*cos(2*pi/24*phase.pc[n])
        amp.b = amp*sin(2*pi/24*phase.pc[n])
        AA= (amp.a -1i*amp.b) * exp(1i*pi/2)
        amp.a = Re(AA)
        amp.b = Im(AA)
        lines(c(0, amp.a), c(0, amp.b), type='l', lwd=1.0, col=bg)
        
        #lines(c(0, tf.a),c(0, tf.b), lwd=1.2, col=bg)
        #text(cof.aa.text[n], cof.bb.text[n], cof.names[n], cex=1.1, col=col, srt=srt, adj=adj)
    }
    legend(-17, -4, legend = c('transcription', 'cytoskeleton', 'chaperone', 'cell cycle', 'DNA repair', 'proteolisis', 'protein transprot', 'others') , cex=1.2, pch=rep(21, 7),
    col='black', pt.bg= bg.cols, text.col='black', pt.cex=1.4, pt.lwd=1.0, border = NA, bty = 'n')
    
    dev.off()

    amps.pc.check = FALSE
    if(amps.pc.check)
    {
        load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v3.Rdata')
        load(file='Tables_DATA/Table_Nuclear_total_mRNA_pol2.Rdata')
        res.c = res.c[, -c(5,6)]
        
        gg = res.m$subunits
        gg = unique(unlist(strsplit(as.character(gg), ',')))
        mm = match(gg, nuclear.names[,3])
        gg = gg[which(!is.na(mm)==TRUE)]
        mm = mm[which(!is.na(mm)==TRUE)]
        amps = nuclear$amp[nuclear.names[mm,1]]
        pvals = nuclear$pval[nuclear.names[mm,1]]
        expression = nuclear.all$mean.mrna[nuclear.names[mm, 1]]
        
        nb.pc = c()
        nb.pc.r = c()
        for(gene in gg)
        {
            nb.pc = c(nb.pc, length(grep(gene, res$subunits)))
            nb.pc.r = c(nb.pc.r, length(grep(gene, res.m$subunits)))
        }
        
        subs = data.frame(gg, nb.pc.r, nb.pc, amps, pvals, expression, stringsAsFactors=FALSE)
        kk = which(subs[,5]<0.05)
        
        nb.c = c()
        amps = c()
        for(n in 1:nrow(res.m))
        {
            gg = res.m$subunits[n]
            gg = unlist(strsplit(as.character(gg), ','))
            mm = match(gg, subs[,1])
            mm = mm[which(!is.na(mm))]
            nb.c = c(nb.c, mean(subs[mm, 3]))
            amps = c(amps, mean(subs[mm, 4]))
        }
        
        plot(subs[,2], subs[,3], cex=0.5)
        abline(0, 1, lwd=2.0, col='red')
        
        plot(subs[kk,4], subs[kk,3], cex=0.5, col='black')
        points(subs[kk, 4], subs[kk, 2], cex=0.5, col='blue')
        
        pdf.name = paste("myplots/FIGURES/low_amps_Rhythmic_Complexes_1.pdf", sep='')
        pdf(pdf.name, width=2.5, height=2.5)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        boxplot(subs[,4] ~ subs[, 3], xlim=c(0, 25), ylim=c(0, 2.0), axes=FALSE, pars = list(boxwex = 0.7, staplewex = 0.4, outwex = 0.4))
        
        axis(1,at=seq(0,25,by=5),cex.axis =1.)
        axis(2,at=seq(0, 2, by=0.5), las=1,cex.axis = 1.)
        #box()
        
        dev.off()
        
        pdf.name = paste("myplots/FIGURES/low_amps_Rhythmic_Complexes_2.pdf", sep='')
        pdf(pdf.name, width=2.0, height=2.0)
        par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        plot(nb.c, amps, cex=0.5, axes=FALSE,xlab=NA, ylab=NA,main=NA)
        axis(1,at=seq(0,50,by=10),cex.axis =1.2)
        axis(2,at=seq(0, 2.5, by=0.5), las=1,cex.axis = 1.2)
        box()
        
        dev.off()
        
    }

    #############
    ###### Examples of Rhythmic protein complexes
    #############
    Plot.PC.examples = TRUE
    if(Plot.PC.examples)
    {
        ##### Examples of rhythmic protein complexes
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
        #Plots.complexes.v2(aa, c(1:nrow(aa)), 'myplots/All_PCs_Corum_subunits_Selected_above_background_v4_rhythmic.pdf')
        #save(res, res.sel, res.m, res.c, file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v2.Rdata' )
        #load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v2.Rdata' )
        load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v3.Rdata')
        res.c = res.c[, -c(5,6)]
        kk = which(res.c$amp.pc<=2.1)
        res.c = res.c[kk,]
        
        aa = res.m
        xx = res.c
        o1 = order(-xx$amp.pc)
        xx = xx[o1,]
        yy = xx
        #indexx = c(1:nrow(res.c))
        pc.names = c('Mediator/asscociated','NCOR/SMRT', 'NCOR-SIN3','NuRD', 'Mt3-Hsp84-Ck', 'HSP90-associated',
                    'CSN/CSA/DDB2', 'APPBP1-UBA3', 'Profilin 1', 'AQP2-force-generator', 'Retromer', 'Endocytic coat',
                    'AMPK', 'CKII', 'DNMT1-RB1-HDAC1-E2F1', 'APC/C')
        mm = match(pc.names, xx$names.pc)
        
        xx = yy[mm, ]
        #pc.names = res.c[,1]
        #pc.cols = c(rep('magenta', 3), rep('darkblue', 5), rep('black', 2))
        
        for(n in 1:nrow(xx))
        {
            name = xx[n, 1]
            jj = which(res.m[,1]==name)
            
            subunits = unique(unlist(strsplit(as.character(res.m$subunits[jj]), ',')))
            index = c()
            for(ss in subunits)
            {
                index = c(index, nuclear.names[which(nuclear.names[,3]==ss) ,1])
            }
            index = unique(index[which(!is.na(index)==TRUE)])
            #kk = unique(c(match(res.m$Complex.id[jj], res.sel$Complex.id), match(res.m$Complex.name[jj], res.sel$Complex.name)))
            #ii = unique(kk[which(!is.na(kk)==TRUE)])
            #if(length(ii)!=length(kk)) print(n)
            #kk = ii
            pc.name = c(as.character(xx[n, 5]), as.character(xx[n, 2]))
            pc.name = unlist(strsplit(as.character(pc.name), ' '))
            pc.name = unlist(strsplit(as.character(pc.name), '/'))
            pc.name = unlist(strsplit(as.character(pc.name), '-'))
            pc.name = paste(pc.name, sep='', collapse='_')
            pc.name = paste(pc.name, '_', signif(xx$amp.pc[n],d=2), sep='', collapse='_')
            
            #genes = unlist(strsplit(as.character(res.sel$subunits.detected[kk]), ','))
            kk = which(nuclear$nb.timepoints[index]>=8)
            index = index[kk]
            #genes = genes[kk]
            
            if(length(index)>1 & xx$amp.pc[n]>0.05)
            {
                pdf.name = paste("myplots/FIGURES/", pc.name, "_Rhythmic_Complexes_examples.pdf", sep='')
                pdf(pdf.name, width=1.8, height=1.6)
                par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
                
                test = as.matrix(nuclear[index,c(1:16)])
                test =apply(test, 1, standadization.nona.impute)
                genes = nuclear$Gene.names[index]
                #test =apply(test, 1, standadize.nona)
                test = t(test)
                pvals = signif(nuclear$pval[index],d=2)
                phases = signif(nuclear$phase[index],d=2)
                amp = signif(nuclear$relamp[index],d=1)
                o1 = order(pvals)
                pvals = pvals[o1]
                test = test[o1,]
                phases = phases[o1]
                index = index[o1]
                genes = genes[o1]
                ylim = range(test, na.rm=TRUE)
                xlim = c(0,48)
                
                cex = 0.8
                plot(1,1,type = 'n', xlab = NA, ylab = NA, xlim = xlim, ylim = ylim, main=paste(xx$names.pc[n]), axes=FALSE, cex.lab=1.0, cex.main=cex)
                abline(h=0, col='darkgray',lwd=1.5, lty=2)
                
                #text.col = c(rep('cornflowerblue', 5),  rep('darkorchid',5))
                #text.col = rep('blue', length(genes))
                text.col = 'black'
                for(i in 1:nrow(test))
                {
                    i.rel = i;
                    type.plot='b';
                    #col='blue'
                    col= 'gray'
                    
                    points(c(0:15)*3, test[i,], type='l', lwd=1.2, col = col , lty = 1)
                    #legend(x = 42,y = (ylim[2]-(i-1)*0.2), legend = paste(genes[i]), text.col=text.col, bty = 'n', cex=0.6)
                }
                bg = 'black'
                if(xx$functions[n]=='Transcriptional regulation') bg = bg.cols[1]
                if(xx$functions[n]=='Cytoskeleton') bg = bg.cols[2]
                if(xx$functions[n]=='Chaperone') bg = bg.cols[3]
                if(xx$functions[n]=='Cell cycle') bg = bg.cols[4]
                if(xx$functions[n]=='DNA repair') bg = bg.cols[5]
                if(xx$functions[n]=='Proteolisis') bg = bg.cols[6]
                if(xx$functions[n]=='Transport') bg = bg.cols[7]
                
                median = apply(test, 2, mean)
                points(c(0:15)*3, median, type='l', col=bg, lwd=2.0)
                axis(1,at=seq(0, 48, by=12),cex.axis =1.)
                axis(2,at = seq(-4, 4, by=1),las=1,cex.axis = 1.)
                #abline(v=aa$phase.p1[kk], col='gray40',lwd=1.2)
                box()
                dev.off()
            }
            
        }
    }
    
}

#################
####### FIGURE 4 (Main + SUpplementary)
#################
FIGURE_4 = TRUE
if(FIGURE_4)
{
    
    ########################
    #### Analysis rhythmic nuclear proteins for rRNA synthesis, processing, translation and assembles
    ########################
    
    #### PolI complex
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
    load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v3.Rdata')
    res.c = res.c[, -c(5,6)]
    kk = which(res.c$amp.pc<=2.1)
    res.c = res.c[kk,]
    
    aa = res.m
    mm = grep('Pol', res.c[,2])
    res.c[mm, ]
    
    mm = 54
    xx = res.c[mm, ]
    #pc.names = res.c[,1]
    #pc.cols = c(rep('magenta', 3), rep('darkblue', 5), rep('black', 2))
    pc.name = 'RNA PolI'
    
    for(n in 1:nrow(xx))
    {
        name = xx[n, 1]
        jj = which(res.m[,1]==name)
        
        subunits = unique(unlist(strsplit(as.character(res.m$subunits[jj]), ',')))
        index = c()
        for(ss in subunits)
        {
            index = c(index, nuclear.names[which(nuclear.names[,3]==ss) ,1])
        }
        index = unique(index[which(!is.na(index)==TRUE)])
        #kk = unique(c(match(res.m$Complex.id[jj], res.sel$Complex.id), match(res.m$Complex.name[jj], res.sel$Complex.name)))
        #ii = unique(kk[which(!is.na(kk)==TRUE)])
        #if(length(ii)!=length(kk)) print(n)
        #kk = ii
        pc.name = c(as.character(xx[n, 5]), as.character(xx[n, 2]))
        pc.name = unlist(strsplit(as.character(pc.name), ' '))
        pc.name = unlist(strsplit(as.character(pc.name), '/'))
        pc.name = unlist(strsplit(as.character(pc.name), '-'))
        pc.name = paste(pc.name, sep='', collapse='_')
        pc.name = paste(pc.name, '_', signif(xx$amp.pc[n],d=2), sep='', collapse='_')
        
        #genes = unlist(strsplit(as.character(res.sel$subunits.detected[kk]), ','))
        kk = which(nuclear$nb.timepoints[index]>=8)
        index = index[kk]
        #genes = genes[kk]
        
        if(length(index)>1 & xx$amp.pc[n]>0.05)
        {
            pdf.name = paste("myplots/FIGURES/", pc.name, "_Rhythmic_Complexes_examples.pdf", sep='')
            pdf(pdf.name, width=2.0, height=1.7)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            test = as.matrix(nuclear[index,c(1:16)])
            test =apply(test, 1, standadization.nona.impute)
            genes = nuclear$Gene.names[index]
            #test =apply(test, 1, standadize.nona)
            test = t(test)
            pvals = signif(nuclear$pval[index],d=2)
            phases = signif(nuclear$phase[index],d=2)
            amp = signif(nuclear$relamp[index],d=1)
            o1 = order(pvals)
            pvals = pvals[o1]
            test = test[o1,]
            phases = phases[o1]
            index = index[o1]
            genes = genes[o1]
            ylim = range(test, na.rm=TRUE)
            xlim = c(0,48)
            plot(1,1,type = 'n', xlab = NA, ylab = NA, xlim = xlim, ylim = ylim, main=NA, axes=FALSE)
            abline(h=0, col='darkgray',lwd=1.5, lty=2)
            
            #text.col = c(rep('cornflowerblue', 5),  rep('darkorchid',5))
            #text.col = rep('blue', length(genes))
            text.col = 'black'
            for(i in 1:nrow(test))
            {
                i.rel = i;
                type.plot='b';
                #col='blue'
                col= 'mediumpurple'
                
                points(c(0:15)*3, test[i,], type='l', lwd=1.2, col = col , lty = 1)
                #legend(x = 42,y = (ylim[2]-(i-1)*0.2), legend = paste(genes[i]), text.col=text.col, bty = 'n', cex=0.6)
            }
            bg = 'black'
            if(xx$functions[n]=='Transcriptional regulation') bg = bg.cols[1]
            if(xx$functions[n]=='Cytoskeleton') bg = bg.cols[2]
            if(xx$functions[n]=='Chaperone') bg = bg.cols[3]
            if(xx$functions[n]=='Cell cycle') bg = bg.cols[4]
            if(xx$functions[n]=='DNA repair') bg = bg.cols[5]
            if(xx$functions[n]=='Proteolisis') bg = bg.cols[6]
            if(xx$functions[n]=='Transport') bg = bg.cols[7]
            
            median = apply(test, 2, mean)
            points(c(0:15)*3, median, type='l', col=bg, lwd=2)
            axis(1,at=seq(0, 48, by=12),cex.axis =1.)
            axis(2,at = seq(-4, 4, by=1),las=1,cex.axis = 1.)
            #axis(1,at=seq(0, 48, by=12),cex.axis =1.2)
            #axis(2,at = seq(-4, 4, by=1),las=1,cex.axis = 1.2)
            #abline(v=aa$phase.p1[kk], col='gray40',lwd=1.2)
            box()
            dev.off()
        }
        
    }

    ###### subgroups of nuclear proteins involved in rRNA synthesis and processing
    for(nn in 1:6)
    {
        eval(parse(text = paste("yy = read.csv('Tables_DATA/rRNA_synthesis_", nn, ".csv', header=TRUE)", sep='')))
        
        #kk = which(!is.na(yy[,20]))
        #yy = yy[kk,]
        xx = as.matrix(yy[, c(1:16)])
        res = t(apply(xx, 1, f24_R2_alt2, t=c(0:15)*3))
        kk = which(res[,1]>=6)
        res = res[kk,]
        xx = xx[kk,]
        
        xx = as.matrix(xx)
        if(ncol(xx)>1)
        {
            xx = t(apply(xx, 1, standadize.nona))
        }else{
            xx = t(xx)
            xx = t(apply(xx, 1, standadize.nona))
        }
        print(length(kk))
        
        pdfname = paste('myplots/FIGURES/rRNA_biogenesis_', nn, '.pdf', sep='')
        
        if(nn==0)
        {
            pdf(pdfname, width=2.2, height=2.2)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            phases = res[,5]
            circular_phase24H_histogram(phases, col='darkgray', cex.axis=0.7, cex.lab=0.01, lwd=1.0)
        }else{
            
            pdf(pdfname, width=2.0, height=1.7)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            lims = range(xx, na.rm=TRUE)
            
            plot(1, 1, type='n', xlim=c(0, 48), ylim=lims, axes=FALSE, xlab=NA, ylab=NA)
            col= 'mediumpurple'
            for(i in 1:nrow(xx))
            {
                
                points(c(0:15)*3, xx[i,], type='l', lwd=1.2, col = col , lty = 1)
                #legend(x = 42,y = (ylim[2]-(i-1)*0.2), legend = paste(genes[i]), text.col=text.col, bty = 'n', cex=0.6)
            }
            bg = 'black'
            
            median = apply(xx, 2, mean, na.rm=TRUE)
            points(c(0:15)*3, median, type='l', col=bg, lwd=2.)
            axis(1,at=seq(0, 48, by=12),cex.axis =1.)
            axis(2,at = seq(-4, 4, by=1),las=1,cex.axis = 1.)
            abline(h=0, col='darkgray',lwd=1.5, lty=2)
            box()
            
        }
        
        dev.off()
    }
    
    ########################
    #### Analysis rhythmic nuclear proteins for DNA damage repair
    ########################
    #### All detected rhythmic proteins
    library(plotrix)
    library("circular")
    make_circ_coord = function(t,x,ttot=24) {
        dt=(t[2]-t[1])*.45
        a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
        h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
        list(angles=a,heights=h)
    }
    circular_phase24H_histogram<-function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5){
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        br=0:24
        h=hist(x, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        radial.plot(co$heights,co$angle,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
    }
    
    phases = xx$phase
    
    pdf('myplots/FIGURES/Phase_distribution_rhythmic_proteins_rRNA_synthesis.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.6,0.2,0.5), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    
    ###### subgroups
    for(nn in 0:8)
    {
        if(nn<2)
        {
            eval(parse(text = paste("yy = read.csv('Tables_DATA/DNA_repair_", nn, ".csv', header=TRUE)", sep='')))
        }else{
            eval(parse(text = paste("yy = read.csv('Tables_DATA/DNA_repair_", nn, ".csv', header=FALSE)", sep='')))
        }
        #kk = which(!is.na(yy[,20]))
        #yy = yy[kk,]
        xx = as.matrix(yy[, c(1:16)])
        res = t(apply(xx, 1, f24_R2_alt2, t=c(0:15)*3))
        kk = which(res[,1]>=6)
        res = res[kk,]
        xx = xx[kk,]
        
        xx = as.matrix(xx)
        if(ncol(xx)>1)
        {
            xx = t(apply(xx, 1, standadize.nona))
        }else{
            xx = t(xx)
            xx = t(apply(xx, 1, standadize.nona))
        }
        print(length(kk))
        
        pdfname = paste('myplots/FIGURES/DNA_damage_repair_', nn, '.pdf', sep='')
        
        if(nn==0)
        {
            pdf(pdfname, width=2.2, height=2.2)
            par(cex = 0.6, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            phases = res[,5]
            circular_phase24H_histogram(phases, col='darkgray', cex.axis=0.7, cex.lab=0.01, lwd=1.0)
        }else{
            
            pdf(pdfname, width=2.0, height=1.7)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            lims = range(xx, na.rm=TRUE)
            
            plot(1, 1, type='n', xlim=c(0, 48), ylim=lims, axes=FALSE, xlab=NA, ylab=NA)
            col= 'cornflowerblue'
            for(i in 1:nrow(xx))
            {
                
                points(c(0:15)*3, xx[i,], type='l', lwd=1.2, col = col , lty = 1)
                #legend(x = 42,y = (ylim[2]-(i-1)*0.2), legend = paste(genes[i]), text.col=text.col, bty = 'n', cex=0.6)
            }
            bg = 'black'
            
            median = apply(xx, 2, mean, na.rm=TRUE)
            points(c(0:15)*3, median, type='l', col=bg, lwd=2.5)
            axis(1,at=seq(0, 48, by=12),cex.axis =1.)
            axis(2,at = seq(-4, 4, by=1),las=1,cex.axis = 1.)
            abline(h=0, col='darkgray',lwd=1.5, lty=2)
            box()

        }
        
        dev.off()
    }

}

#################
####### FIGURE 5 (Main + SUpplementary)
#################
FIGURE_5 = TRUE
if(FIGURE_5)
{
    load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Rdata/Phos_Nuclear_Mapping_all_v2.Rdata')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    
    #write.table(phospho.nuclear, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Table_SXX_nuclear_phospho_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v3.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
    nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
    nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')
    #load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Table_Nuclear_mRNA_Total_Nascent.Rdata')
    load(file='Tables_DATA/Table_Nuclear_total_mRNA_pol2.Rdata')
    
    Table.Sup = FALSE
    if(Table.Sup)
    {
        xx = phospho.nuclear
        yy = nuclear.all
        mm = match(xx$mapping.gene, yy$Gene.names)
        xx = cbind(xx, yy[mm, c(76:95)])
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_SXX_phospho_nuclear_proteins_mRNAs_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        xx = phospho.nuclear
        kk = c(42, 40, 32, 34, 36, 37:39, 41, 1:16, 43:49, 51:73)
        xx = xx[, kk]
        colnames(xx)[c(10:25)] = paste('ZT', c(0:15)*3, '.phospho', sep='')
        colnames(xx)[c(26:32)] = paste(colnames(xx)[c(26:32)], '.phospho', sep='')
        #xx = nuclear[, c(20, 19, 17, 18, 1:16, 34:37, 22:28)]
        #colnames(xx)[25:31] = paste(colnames(xx)[25:31], '.WT', sep='')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S5_phospho_nuclear_proteins_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        load(file='Tables_DATA/Kinases_Candidates_list_curated_targets.Rdata')
        
        xx = kinase.candidates[,-1]
        xx = xx[, c(1, 2, 4, 3, 12)]
        
        xx = data.frame(rownames(xx), xx, stringsAsFactors=FALSE)
        colnames(xx)[1:3] = c('kinase.motif', 'phase.LM', 'ampl.LM')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S5_kinases_inferred.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    }

    
    ########
    #### Global analysis of phosphorylated sites
    ########
    #########
    ##### Global comparison between phospho and nuclear proteins
    aa = phospho.nuclear
    
    kk = which(phospho$nb.timepoints>=8)
    gene.names = phospho$gene[kk]
    gene.names = gene.names[which(gene.names!=' ')]
    length(unique(gene.names))
    
    ### overlap between phospho and total
    kk = which(phospho$nb.timepoints>=8 & phospho$gene!='')
    jj = which(nuclear$nb.timepoints>=8 & nuclear$Gene.names!='')
    
    require('VennDiagram')
    cols = c( "lightskyblue1", "lightgreen")
    venn.diagram(list(phophos = phospho$gene[kk], nuclear = nuclear$Gene.names[jj]),
    height = 2.0, width = 2.0, units = "in",
    fill = cols, cat.col='black', resolution=500,
    alpha = c(0.6, 0.6), cex = 1.0, cat.fontface = 2,lty=1,sub.cex=0.8,
    mar = rep(0., 4)+0.1, cat.dist=c(0.1, 0.08), ext.dist=0.05,
    filename = "myplots/FIGURES/Phospho_nuclear_overlap.png");
    
    ### rhythms in phospho
    qq = c(0:100)/100
    nb.rhythmic = c()
    for(n in 1:length(qq))
    {
        cutoff.qq = qq[n]
        nb.rhythmic = c(nb.rhythmic,length(which(phospho$qv<cutoff.qq)))
    }
    
    pdf('myplots/FIGURES/Nb_rhythmic_phosphos_proteins_FDR.pdf', width=2.0, height=2.0)
    par(cex = 0.7, las = 0, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(qq, nb.rhythmic+0.01, type='l', lty=1, lwd=1.5, main='', log='y',ylim=c(1,5000), xlab=NA, ylab=NA, col=c('black'), axes=FALSE)
    #points(qq, nb.rhythmic.robles+0.01, type='l', lty=1, lwd=2, col='green')
    abline(v=0.05,col='darkgreen',lwd=2.0,lty=2)
    #abline(v=0.2,col='darkgreen',lwd=2,lty=2)
    #abline(v=0.1,col='darkgreen',lwd=2,lty=2)
    #legend('topright', legend = c('Mauvoisin','Robles'), lty=c(1,1), cex=0.7,col = c('blue', 'green'), border = NA, bty = 'n')
    axis(1,at=c(0.05,0.2, 0.5, 1),cex.axis =1.2)
    axis(2,at = c(1, 10, 100, 1000),las=1,cex.axis = 1.2)
    #abline(v=aa$phase.p1[kk], col='gray40',lwd=1.2)
    box()

    dev.off()
    
    type= c('S', 'T','Y')
    jj = which(phospho$nb.timepoints>=8)
    coverages = c(length(which(phospho$Amino.acid[jj]=='S')), length(which(phospho$Amino.acid[jj]=='T')),
                    length(which(phospho$Amino.acid[jj]=='Y')))
    #kk = which(phospho$qv<0.05 & phospho$nb.timepoints>=8)
    rhythm = c(length(which(phospho$Amino.acid[jj]=='S' & phospho$qv[jj]<0.05)), length(which(phospho$Amino.acid[jj]=='T' & phospho$qv[jj]<0.05)), length(which(phospho$Amino.acid[jj]=='Y' & phospho$qv[jj]<0.05)))
    rhythm = rhythm/coverages
    coverages = coverages/length(jj)
    counts <- rbind(coverages, rhythm)
    stat = t(counts)
    coverages = stat[,1]*(stat[,2])
    rhythm = stat[,1]*(1-stat[,2])
    counts <- rbind(coverages, rhythm)*100

    pdf.name = paste("myplots/FIGURES/Phospho_sites_rhythmic_percentage.pdf", sep='')
    pdf(pdf.name, width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(5,3,2,0.8)+0.1, tcl = -0.3)
    
    barplot(counts, main=NA, col=c("gray40", "gray"), ylim=c(0, 110), space=0.5,border=FALSE,
    legend.text=c('', ''), las=2, cex.lab=1., args.legend=list(x=3, y=120, bty='n', cex=1, border=FALSE), axes=FALSE)
    #abline(h=0.15, col='red', lwd=1.5)
    axis(2,at = seq(0, 100, by=20),las=1,cex.axis = 1.)
    dev.off()
    
    
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
    
    #### Rhythmic phospho peptides: phases and amplitudes
    jj = which(phospho$qv<0.05)
    phases = phospho$phase[jj]
    amps = phospho$amp[jj]
    
    pdf('myplots/FIGURES/Phase_distribution_rhythmic_phospho_peptides_detected.pdf', width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col='darkgreen', cex.axis=0.7, cex.lab=0.01, lwd=0.6)
    
    dev.off()
    
    pdf('myplots/FIGURES/Amplitudes_distribution_rhythmic_phospho_peptides_detected.pdf', width=2.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=c(0:12)/2;
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=1.0)
    axis(1)
    axis(2)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=1.0)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=1.0)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=1.0)
    lines(h$breaks, rep(0,length(h$breaks)), type='S')
    invisible(h)
    
    dev.off()
    
    #### Compare WT and CRY DKO and annotate TFs and cofactors
    ##### WT vs. Cry1/2 KO
    bb = phospho[which(phospho$qv<0.05), ]
    res.ko = c()
    for(n in 1:nrow(bb))
    {
        #index = n
        data.wt = as.numeric(bb[n, c(1:16)]);
        res.ko = rbind(res.ko, chow.test(c(data.wt, as.numeric(bb[n, grep('Cry.KO', colnames(bb))]))))
    }
    bb$pval.cryko = res.ko[,1]
    bb$nb.cryko = res.ko[,2]
    bb$cor.cryko = res.ko[,3]
    
    length(which(bb$cor.cryko>0.5))/nrow(bb)
    length(which(bb$pval.cryko>0.05))/nrow(bb)
    length(which(bb$pval.cryko>0.05 & bb$cor.cryko>0.5))/nrow(bb)
    
    cc = nuclear[which(nuclear$qv<0.05), ]
    res.ko = c()
    for(n in 1:nrow(cc))
    {
        data.wt = as.numeric(cc[n, c(1:16)]);
        res.ko = rbind(res.ko, chow.test(c(data.wt, as.numeric(cc[n, grep('Cry.KO', colnames(cc))]))))
    }
    cc$pval.cryko = res.ko[,1]
    cc$nb.cryko = res.ko[,2]
    cc$cor.cryko = res.ko[,3]
    
    length(which(cc$cor.cryko>0.5))/nrow(cc)
    length(which(cc$pval.cryko>0.05))/nrow(cc)
    length(which(cc$pval.cryko>0.05 & cc$cor.cryko>0.5))/nrow(cc)
    
    
    ##### Phospho vs. Nuclear
    kk = which(aa$nb.timepoints>=6 & aa$nb.timepoints.nuclear>=6)
    aa = aa[kk, ]
    names = aa$gene
    
    res.ko = c()
    for(n in 1:nrow(aa))
    {
        #index = n
        data.wt = as.numeric(aa[n, c(1:16)]);
        res.ko = rbind(res.ko, chow.test(c(data.wt, as.numeric(aa[n, grep('Cry.KO', colnames(bb))]))))
    }
    aa$pval.cryko = res.ko[,1]
    aa$nb.cryko = res.ko[,2]
    aa$cor.cryko = res.ko[,3]
    
    tfs.cofactors = TRUE
    if(tfs.cofactors)
    {
        mapping.all = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Manually_curation_TFs_pools.txt', header=TRUE, sep='\t')
        tfs.all = mapping.all$TFs.pool
        mcors = c('Ncor1', 'Ep300', 'Hdac1', 'Hdac2', 'Hdac5', 'Sirt6', 'Crebbp', 'Ncoa1', 'Ncoa2', 'Ncoa3', 'Ncor1', 'Ncor2', 'Nrip1', 'Kat2a', 'Kat2b', 'Hdac1', 'Hdac2',
        'Hdac3', 'Hdac4', 'Hdac5', 'Hdac6', 'Hdac7', 'Hdac8', 'Hdac9', 'Crtc1', 'Crtc2', 'Crtc3', 'Rb1', 'Med1', 'Ppargc1a', 'Ppargc1b', 'Setd2', 'Setdb1', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2')
        tfs.all = tfs.all[which(is.na(match(tfs.all, mcors))==TRUE)]
        cofactor = read.delim('Annotations/Cofactors_goterms.txt',header=FALSE, sep='\t', as.is = c(3));
        cofactor = unique(c(cofactor[,3], mcors))
        
        index = c()
        phospho.rhy = c()
        for(n in 1:length(names))
        {
            gene = unlist(strsplit(as.character(names[n]), ';'))
            index = c(index, rep(n, length(gene)))
            phospho.rhy = c(phospho.rhy, gene)
        }
        index = data.frame(index, phospho.rhy, stringsAsFactors=FALSE)
        
        tfs.cors = rep(NA, nrow(index))
        
        jj = match(index[,2], tfs.all)
        mm = match(index[,2], cofactor)
        jj = which(!is.na(jj)==TRUE)
        mm = which(!is.na(mm)==TRUE)
        
        index[intersect(jj,mm),2]
        
        tfs.cors[jj] = 'tfs'
        tfs.cors[mm] = 'cors'
        index = data.frame(index, tfs.cors, stringsAsFactors=FALSE)
    }
    
    ##########
    ##### Two groups of genes: nuclear flat, phospho rhythmic; nuclear rhythmic- phospho rhythmic
    ##########
    
    #########################
    ### Class I:: Rhythmic nuclear and Rhythmic phospho
    
    ### Heatmap
    fdr.cutoff = 0.05
    jj = which(aa$qv<fdr.cutoff & aa$qv.nuclear<fdr.cutoff)
    
    x=aa[jj, c(51:66)]
    y=aa[jj, c(1:16)]
    
    #x = t(apply(x, 1, centering.nona))
    #y = t(apply(y, 1, centering.nona))
    
    o=order(aa$phase[jj])
    x=x[o,]
    x = t(apply(x, 1, standadize.nona))
    y=y[o,]
    y = t(apply(y, 1, standadize.nona))
    #x = cbind(x,rep(NA,nrow(x)))
    xy=cbind(x,y)
    
    library("gplots")
    xy = as.matrix(xy)
    x = as.matrix(x)
    y = as.matrix(y)
    
    pdf('myplots/FIGURES/heatmap_TFs_Cors_phospho_nuclear_R_R_nuclear.pdf', width=1.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    #heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0, ColSideColors=c(rep('gray70',4),rep('black',4),rep('gray70',4),rep('black',4)), margins = c(0.1, 1.7),lhei = c(0.1,4.0),lwid=c(0.05,2.5),
    #		  labRow=NA,cexRow=0.5,offsetRow=0.1, xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
    
    heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0, margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
    labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
    
    dev.off()
    
    pdf('myplots/FIGURES/heatmap_TFs_Cors_phospho_nuclear_R_R_phospho.pdf', width=1.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0, margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
    labRow=NA, cexRow=0.55,offsetRow=0.1, xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
    
    dev.off()
    
    ### phases and amplitudes
    fdr.cutoff = 0.05
    jj = which(aa$qv<fdr.cutoff & aa$qv.nuclear<fdr.cutoff)
    phases = aa$phase[jj]
    amps = aa$amp[jj]
    
    pdf('myplots/FIGURES/Phase_distribution_rhythmic_phospho_Class_I_RR.pdf', width=1.7, height=1.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col='gray', cex.axis=0.6, cex.lab=0.01, lwd=0.6)
    
    dev.off()
    
    pdf('myplots/FIGURES/Amplitudes_distribution_rhythmic_phospho_Class_I_RR.pdf', width=1.6, height=1.6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=c(0:12)/2;
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=1.0)
    axis(1)
    axis(2)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=1.0)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=1.0)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=1.0)
    lines(h$breaks, rep(0,length(h$breaks)), type='S')
    invisible(h)
    
    dev.off()

    
    ### Examples of Class I and II
    index.aa = c()
    
    for(n in 1:nrow(aa))
    {
        index.aa = c(index.aa, paste(c(aa$gene[n], as.character(aa$Amino.acid[n]), aa$Position[n]), sep='', collapse='_'))
    }
    
    #fdr.cutoff = 0.05
    examples = c('Rorc_T_119', 'Per2_S_697', 'Dbp_S_86', 'Hdac3_S_424', 'Hnf1b_S_80', 'Elf1_S_187', 'Foxk1_S_431', 'Hbp1_S_46',
                'Nfia_S_280', 'Foxp2_S_328', 'Cic_S_1080', 'Chd4_S_1517', 'Ncbp1_S_22')
    jj = match(examples, index.aa)
    
    
    for(n in jj)
    {
        gene = aa$gene[n]
        
        print(gene)
        pdf.name = paste("myplots/FIGURES/Phospho_examples_Nuclear_centered_GroupS_RR_RF_",toupper(index.aa[n]),".pdf", sep = "")
        pdf(pdf.name, width=1.8, height=1.6)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        y0 = as.numeric(aa[n, c(1:16)])
        y1 = as.numeric(aa[n, c(51:66)])
        y2 = as.numeric(aa[n, grep('Cry.KO', colnames(aa))])
        y2 = y2-mean(y0, na.rm=TRUE)
        y0 = y0-mean(y0, na.rm=TRUE)
        y1 = y1 - mean(y1, na.rm=TRUE)
        
        #averg0 = mean.err(y0)[1,]
        #err0 = mean.err(y0)[2,]
        #averg1 = mean.err(y1)[1,]
        #err1 = mean.err(y1)[2,]
        
        time = c(0:15)*3
        lims = range(c(y0, y1), na.rm=TRUE)
        
        xlim = c(0, 48)
        col = 'darkgreen'
        pch = 16
        lwd = 1.2
        cex = 1.0
        plot(time, y0, type='l', ylim=lims, xlim=xlim, col=col, lwd=lwd,  pch=pch, lty=1, axes=FALSE, xlab=NA, ylab=NA,  main=paste(toupper(gene), ', ', aa$Amino.acid[n], aa$Position[n], sep=''), cex.main=0.8)
        #arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
        points(time, y0, type='p', col=col, pch=pch, cex=0.7)
        col='indianred1'
        #points(c(0,6,12, 18), y2, type='b', col=col, lwd=lwd, lty=2, pch=pch)
        
        col='darkblue'
        lwd = 1.2
        #pch = 15
        points(time, y1, type='l', col=col, lwd=lwd, lty=1, pch=pch)
        #arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3,  lty=1, col=col, pch=pch,lwd=lwd)
        points(time, y1, type='p', col=col, pch=pch, cex=0.7)
        
        axis(1,at=c(0:4)*12,cex.axis =cex)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        #lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
        
        diff = lims[2]-lims[1]
        if(diff>=3. & diff<6.0) axis(2,at = seq(-5, 5, by=1),las=1,cex.axis = cex)
        if(diff>=1.0 & diff<3.) axis(2,at = seq(-5, 5, by=0.5),las=1,cex.axis = cex)
        if(diff<1) axis(2,at = seq(-5, 5, by=0.2),las=1,cex.axis = 1.)
        
        box()
        
        #axis(2, at= seq(2000, 5000, by=1000), las=0,cex.axis = 0.7)
        #box()
        #legend(30, 4600, c('all','mononu', 'binu'), lty=c(1,2,4),col=c('blue', 'darkgray', 'darkgray'), cex=0.7, pt.lwd=1.5, pt.cex=1.0, bty='n')
        
        abline(h=0,lty=2,lwd=2.0, col="darkgray")
        #if(gene=='Rorc')
        #legend(0,(lims[2]+0.5), cex=1.0, c('phospho','nuclear' ), pch=c(16,16),lty=c(1,1),col=c("darkgreen", "darkblue"), bty='n')
        
        dev.off()
        
    }
    
    ###### Compare phases and relative amplitudes between phospho and nuclear for Class I
    TFs.check = FALSE
    if(TFs.check)
    {
        kk = which(xx$tfs.cors=='tfs')
        yy = xx[kk,]
        kk = yy[,1]
        
        yy = cbind(yy, aa$qv[kk], aa$qv.nuclear[kk], aa$phase[kk], aa$phase.nuclear[kk], aa$relamp[kk], aa$relamp.nuclear[kk], aa$gene[kk])
        colnames(yy) = c('index', 'gene', 'tfs.cors', 'qv', 'qv.nuclear', 'phase', 'phase.nuclear', 'relamp', 'relamp.nuclear', 'gene.mapping')
        
        kk = which(yy$qv.nuclear<0.05 & yy$qv>=0.05)
        unique(yy$gene[kk])
        kk = which(yy$qv.nuclear<0.05 & yy$qv<0.05)
        unique(yy$gene[kk])
        
        kk = which((yy$qv.nuclear>=0.05 | is.na(yy$qv.nuclear)==TRUE) & yy$qv<0.05)
        unique(yy$gene[kk])
        
    }
    
    kk = which(index[,2] == 'Nr3c1')
    index[kk, ]
    index[kk, 3] = 'tfs'
    
    cutoff = 0.05
    fdr.cutoff = 0.05
    jj = which(aa$qv<fdr.cutoff & aa$qv.nuclear<fdr.cutoff)
    bb = aa[jj, ]
    
    pdf(paste('myplots/FIGURES/Phases_phospho_nuclear_FDR_', cutoff, '.pdf', sep=''), width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    plot(1, 1, type='n', xlim=c(0, 24), ylim=c(0, 24), xlab=NA, ylab=NA, axes=FALSE)
    abline(0,1, lwd=1.2, col='gray')
    abline(24,1, lwd=1.2, col='gray')
    axis(1,at=seq(0, 24, by=6),cex.axis =1.0)
    axis(2,at = seq(0, 24, by=6),las=1,cex.axis = 1.0)
    box()
    
    for(n in 1:nrow(bb))
    {
        test = which(index[,1]==jj[n])
        test = index[test, 3]
        test = unique(test)
        test = test[which(!is.na(test)==TRUE)]
        
        if(length(test)==0) bg = 'white'
        if(length(test)==1)
        {
            if(test=='tfs') {
                bg='magenta';
                print(bb$gene[n]);
            }
            if(test=='cors') {
                bg = 'cornflowerblue';
                #print(bb$gene[n]);
            }
        }
        if(length(test)>1) print(bb$gene[n])
        #if(bb$tfs.cors[n]=='cors') bg='darkgreen'
        #if(bb$tfs.cors[n]=='tfs')
        
        points(bb$phase.nuclear[n], bb$phase[n], pch=21, col='black', lwd=0.7, cex=0.8, bg=bg)
        #text(bb$phase.nuclear[jj], bb$phase[jj], bb$gene.mapping[jj], cex=0.4,col='darkblue', pos=pos, offset=0.2)
        
    }
    legend(2, 24, legend = c('TFs','Cors', 'others') , cex=0.8, pch=c(21, 21), col='black', pt.cex=1.0, pt.bg= c("magenta", "cornflowerblue", "white"), border = NA, bty = 'n')
    #legend(0, 24, c('TFs','Cors'), pch=c(21,21), col='black', bg=c("magenta", "cornflowerblue"), bty='n')
    
    dev.off()
    
    pdf(paste('myplots/FIGURES/Relamp_phospho_nuclear_FDR_', cutoff, '.pdf', sep=''), width=2.0, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    plot(1, 1, type='n', xlim=c(0, 1), ylim=c(0, 1), xlab=NA, ylab=NA, axes=FALSE)
    abline(0,1, lwd=1.2, col='gray')
    #abline(24,1, lwd=1.2, col='gray')
    axis(1,at=seq(0, 1, by=0.2),cex.axis =1.0)
    axis(2,at = seq(0, 1, by=0.2),las=1,cex.axis = 1.0)
    box()
    
    for(n in 1:nrow(bb))
    {
        test = which(index[,1]==jj[n])
        test = index[test, 3]
        test = unique(test)
        test = test[which(!is.na(test)==TRUE)]
        
        if(length(test)==0) bg = 'white'
        if(length(test)==1)
        {
            if(test=='tfs') {
                bg='magenta';
                print(bb$gene[n]);
            }
            if(test=='cors') {
                bg = 'cornflowerblue';
            }
        }
        if(length(test)>1) print(bb$gene[n])
        
        points(bb$relamp.nuclear[n], bb$relamp[n], pch=21, col='black', lwd=0.7, cex=0.8, bg=bg)
        #text(bb$phase.nuclear[jj], bb$phase[jj], bb$gene.mapping[jj], cex=0.4,col='darkblue', pos=pos, offset=0.2)
        
    }
    #legend(0.2, 1, legend = c('TFs','Cors') , cex=1.0, pch=c(21, 21), col='black', pt.cex=1.0, pt.bg= c("magenta", "cornflowerblue"), border = NA, bty = 'n')
    #legend(0, 24, c('TFs','Cors'), pch=c(21,21), col='black', bg=c("magenta", "cornflowerblue"), bty='n')
    
    dev.off()
    
    #######
    ### Class II flat or non-detected nuclear and rhythmic phospho
    #######
    #### heatmap
    fdr.cutoff = 0.05
    jj = which(aa$qv<fdr.cutoff & aa$qv.nuclear>=fdr.cutoff)
    
    x=aa[jj, c(51:66)]
    y=aa[jj, c(1:16)]
    
    x = t(apply(x, 1, standadize.nona))
    y = t(apply(y, 1, standadize.nona))
    
    o=order(aa$phase[jj])
    x=x[o,]
    y=y[o,]
    #x = cbind(x,rep(NA,nrow(x)))
    xy=cbind(x,y)
    
    library("gplots")
    xy = as.matrix(xy)
    x = as.matrix(x)
    y = as.matrix(y)
    
    pdf('myplots/FIGURES/heatmap_TFs_Cors_phospho_nuclear_R_F_nuclear.pdf', width=1.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    #heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, scale='row', labCol=NA, keysize=1.0, ColSideColors=c(rep('gray70',4),rep('black',4),rep('gray70',4),rep('black',4)), margins = c(0.1, 1.7),lhei = c(0.1,4.0),lwid=c(0.05,2.5),
    #		  labRow=NA,cexRow=0.5,offsetRow=0.1, xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50')
    
    heatmap.2(x,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0, margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
    labRow=NA,cexRow=0.5,xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
    
    dev.off()
    
    pdf('myplots/FIGURES/heatmap_TFs_Cors_phospho_nuclear_R_F_phospho.pdf', width=1.0, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    heatmap.2(y,col=redgreen(75),Rowv = NA, Colv=NA, na.rm = TRUE, labCol=NA, keysize=1.0,margins = c(0.1, 0.1),lhei = c(0.05,4.0),lwid=c(0.05,3.0),
    labRow=NA, cexRow=0.55,offsetRow=0.1, xlab=NA, key=FALSE, density.info="none", symkey=FALSE,trace="none",dendrogram="none",na.color='gray50',colsep=c(24),sepwidth=c(0.5,1), sepcolor="darkblue")
    
    dev.off()
    
    #### Phase distribution of phospho in class II
    library(plotrix)
    library("circular")
    make_circ_coord = function(t,x,ttot=24) {
        dt=(t[2]-t[1])*.45
        a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
        h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
        list(angles=a,heights=h)
    }
    circular_phase24H_histogram<-function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5){
        #color_DHS = rgb(0.6,0,0.2)
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(x, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        radial.plot(co$heights,co$angle,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
    }
    
    fdr.cutoff = 0.05
    jj = which(aa$pval<fdr.cutoff & aa$pval.nuclear>=fdr.cutoff)
    phases = aa$phase[jj]
    amps = aa$amp[jj]
    
    pdf('myplots/FIGURES/Phase_distribution_phospho_Class_II.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col='gray50', cex.axis=0.7, cex.lab=0.01, lwd=1.0)
    
    dev.off()
    
    pdf('myplots/FIGURES/Amplitudes_distribution_rhythmic_phospho_Class_II.pdf', width=2., height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    breaks=c(0:6)/2;
    h = hist(amps, breaks=breaks, plot=FALSE)
    y = h$counts
    y[which(y<=0)] = 0.5;
    plot(h$breaks, c(NA,y), type='S', ylim=c(1, max(y)), main=NA, xlab=NA, ylab=NA, axes=FALSE, log='y', lwd=1.0)
    axis(1, at=c(0, 1, 2, 3), cex=1.0)
    axis(2, cex=1.0)
    lines(h$breaks, c(h$counts,NA), type='s', lwd=1.0)
    lines(h$breaks, c(NA,h$counts), type='h', lwd=1.0)
    lines(h$breaks, c(h$counts,NA), type='h',lwd=1.0)
    lines(h$breaks, rep(0,length(h$breaks)), type='S')
    invisible(h)
    
    dev.off()
    
    ##########################
    ###### Activities of kinase motifs and compared with nuclear or phospho proteins
    ##########################
    load(file='Rdata/Phos_Nuclear_Mapping_all_v2.Rdata')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    
    table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v3.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
    nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
    nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')
    load(file='Rdata/Phos_all_with_seq_names_v2.Rdata')
    
    load(file='Rdata/kinase_motifs_curated_occurrence_matrix_v2.Rdata')
    
    load(file='Rdata/Kinase_candidates_elastic_net_alpha_0.02_curated.Rdata')
    load(file='Rdata/Kinase_candidates_phase_enrichment_p_0.05_curated.Rdata')
    load(file='Rdata/kinases_phosphatases_keep.Rdata')
    
    load(file='Tables_DATA/Kinases_Candidates_list_curated_targets.Rdata')
    
    mat.oc = mat.ocm
    kmotifs = kmotifs
    kinase.candidates = kinase.candidates
    infer = kinase.candidates[,-1]
    
    
    ##### Plot phases of inferred Motifs
    phase.m = infer[,1]
    amp.m = infer[,2]
    o1 = order(phase.m)
    phase.m = phase.m[o1]
    amp.m = amp.m[o1]
    motif.names = rownames(infer)[o1]
    infer = infer[o1,]
    
    yy = phase.m
    for(n in 1:(length(yy)-1))
    {
        diff = (yy[n+1]-yy[n]);
        if(diff>=24) diff = diff-24;
        if(diff<0.2) yy[n+1] = yy[n]+0.2;
        if(yy[n+1]>=24) yy[n+1] = yy[n+1]-24;
    }
    phase.m = yy
    
    motif.amp = 2.0;
    amp = motif.amp
    motif.a = amp*cos(2*pi/24*phase.m)
    motif.b = amp*sin(2*pi/24*phase.m)
    CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
    motif.aa = Re(CC)
    motif.bb = Im(CC)
    
    amp = motif.amp + 0.08
    motif.a = amp*cos(2*pi/24*phase.m)
    motif.b = amp*sin(2*pi/24*phase.m)
    CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
    motif.txt1.aa = Re(CC)
    motif.txt1.bb = Im(CC)
    
    #amp = motif.amp-1.
    #motif.a = amp*cos(2*pi/24*infer$phase.nuclear)
    #motif.b = amp*sin(2*pi/24*infer$phase.nuclear)
    #CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
    #motif.txt2.aa = Re(CC)
    #motif.txt2.bb = Im(CC)
    
    rmm=max(abs(c(motif.aa,motif.bb)))+2.0
    rr=c(-rmm,rmm)
    xlim = rr
    ylim = rr
    
    alpha = 0.02
    
    pdf(paste('myplots/FIGURES/Kinase_Motif_activities_inferred_Elastic_net_curated_alpha_', alpha, '.pdf', sep=''), width=6, height=6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1,1,0.8,0.8)+0.1, tcl = -0.3)
    
    plot(motif.aa, motif.bb, main=NA, type='n', xlim=xlim, ylim=ylim, axes=F, xlab='', ylab='', lwd=2, pch=23, col='black',bg='green',cex=2.0)
    
    #abline(v=0,h=0, col='darkgray',lwd=2.0)
    phi=seq(0,2*pi,len=1000)
    #lines(rm*cos(phi), rm*sin(phi), col='darkgray', lwd=2)
    motif.cycle = motif.amp-0.05
    lines(motif.cycle*cos(phi), motif.cycle*sin(phi), col='darkgray', lwd=2.0)
    #lines((motif.amp-0.1)*cos(phi), (motif.amp-0.1)*sin(phi), col='darkgray', lwd=2)
    lines(motif.cycle*cos(phi)/5*4, motif.cycle*sin(phi)/5*4, col='darkgray', lwd=1.0, lty=2)
    lines(motif.cycle*cos(phi)/5*3, motif.cycle*sin(phi)/5*3, col='darkgray', lwd=1.0, lty=2)
    lines(motif.cycle*cos(phi)/5*2, motif.cycle*sin(phi)/5*2, col='darkgray', lwd=1.0, lty=2)
    lines(motif.cycle*cos(phi)/5*1, motif.cycle*sin(phi)/5*1, col='darkgray', lwd=1.0, lty=2)
    
    arrows((motif.amp-0.1), 0, -(motif.amp-0.1), 0,  col='darkgray', lty=2, lwd=1.0, length=0.0);
    arrows(0, (motif.amp-0.1), 0, -(motif.amp-0.1), col='darkgray', lty=2, lwd=1.0, length=0.0);
    
    rainbow = rainbow(length(motif.names),s = 0.85, v = 0.85)
    for(n in 1:length(motif.names))
    {
        bg = 'purple'
        #if(!is.na(infer$phase.MEA[n])) bg = 'purple4'
        points(motif.aa[n], motif.bb[n], pch=21, cex=1.2, col='black',bg=bg)
        
        if(phase.m[n]<=6) srt = 90-phase.m[n]/24*360;
        if(phase.m[n]>6 & phase.m[n]<=12) srt = 90-phase.m[n]/24*360;
        if(phase.m[n]>12 & phase.m[n]<=18) srt = 270-phase.m[n]/24*360;
        if(phase.m[n]>=18) srt = 270-phase.m[n]/24*360;
        cex = 0.7
        col = 'black'
        
        if(phase.m[n]<=12)
        {
            adj = 0;
            #text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=cex, col='darkblue', srt=srt, adj=adj)
            #adj = 1;
            test = unlist(strsplit(as.character(motif.names[n]), '_'))
            test = paste(test[-length(test)], sep='', collapse='_')

            text(motif.txt1.aa[n], motif.txt1.bb[n], toupper(test), cex=cex, col=col, srt=srt, adj=adj)
        }else{
            adj = 1;
            #text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=cex, col='darkblue', srt=srt, adj=adj)
            test = unlist(strsplit(as.character(motif.names[n]), '_'))
            test = paste(test[-length(test)], sep='', collapse='_')
            
            text(motif.txt1.aa[n], motif.txt1.bb[n], toupper(test), cex=cex, col=col, srt=srt, adj=adj)
        }
        
        Add.amps = TRUE
        if(Add.amps)
        {
            ### add amplitudes of motifs (average amplitudes per motif)
            #amp = ((6+log10(amp.m[n]))/5)*motif.cycle
            amp = (amp.m[n]/0.05)*motif.cycle
            amp.a = amp*cos(2*pi/24*phase.m[n])
            amp.b = amp*sin(2*pi/24*phase.m[n])
            AA= (amp.a -1i*amp.b) * exp(1i*pi/2)
            amp.a = Re(AA)
            amp.b = Im(AA)
            lines(c(0, amp.a), c(0, amp.b), type='l', lwd=1.0, col=bg)
        }
        
    }
    #legend(-3, -1.2, legend = c(' ', ' '), cex=1.3, pt.cex=1., col='black', pt.bg = c('chartreuse4', 'darkolivegreen'), pch=21, bty = 'n' )
    
    #}
    dev.off()
    
    #### Validation some examples by comparing temporal activities with nuclear protein levels
    standardization.nona = function(x)
    {
        xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
        return(xx)
    }

    jj = which(infer$pval.nuclear<0.05)
    infer[jj,]
    examples = rownames(infer)[jj]
    
    for(motif in examples)
    {
        
        print(motif)
        pp = infer$phase[which(rownames(infer)==motif)]
        amp = infer$ampl[which(rownames(infer)==motif)]
        mm = as.numeric(infer$index.nuclear[which(rownames(infer)==motif)])
        gene = nuclear$Gene.names[mm]
        print(gene)
        
        pdf.name = paste("myplots/FIGURES/Inferred_Kinases_vs_nuclear_",gene,".pdf", sep = "")
        pdf(pdf.name, width=2.0, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        y0 = standardization.nona(as.numeric(nuclear[mm, c(1:16)]))
        time = seq(0, 48, by=0.5)
        y1 = amp*cos(2*pi/24*(time-pp))
        y1 = standardization.nona(y1)
        
        lims = range(c(y0, y1), na.rm=TRUE)
        
        col = 'purple'
        #if(motif=="Mapk14_ST"| motif=="Csnk1a1_Limk2_Mapkapk3_ST") col='purple4';
        ttest = motif;
        ttest = unlist(strsplit(as.character(ttest), '_'))
        ttest = paste(ttest[-length(ttest)], sep='', collapse='_')
        
        plot(c(0:15)*3, y0, type='l', col="darkblue",lwd=1.2, main=toupper(gene), cex.main=0.8, ylim=lims, ylab=NA, xlab=NA, xlim=c(0,48), axes=FALSE)
        points(c(0:15)*3, y0, type='p', pch=16, col='darkblue', cex=0.8)
        points(time, y1, type='l', lwd=1.5, col=col)
        
        axis(1,at=seq(0, 48, by=12),cex.axis =1.0)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = signif(lims, d=1)
        if(gene=='Hdac3') {
            axis(2,at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),las=1,cex.axis = 1.0)
        }else{
            by = signif((lims[2]-lims[1])/5,d=0)
            axis(2,at = seq(lims[1], lims[2], by=by),las=1,cex.axis = 1.0)
        }
        box()
        abline(h=0,lty=2,lwd=1.5, col="darkgray")
        #if(gene=='Rorc')
        #legend(0,(lims[2]+0.5),c('phospho','nuclear'),pch=c(16,16),lty=c(1,1),col=c("darkgreen", "darkblue"), bty='n')
        
        dev.off()
    }
    
    #######
    #### plot targets of CKD1, CKD4, CDK6 and Wee1
    #######
    rsum = apply(mat.oc, 1, sum)
    
    mtfs = c("Csnk1a1_Limk2_Mapkapk3_ST", "Gsk3a_ST", "Mapk14_ST", "Cdk1_ST", "Cdk3_4_5_6_ST", "Csnk1g3_Rps6ka4_Wee1_ST")
    kk = match(mtfs, rownames(infer))
    
    ref = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Database_PhosphoNetworks/rawKSI.csv',sep='\t', header=TRUE)
    
    for(n.index in 1:length(kk))
    {
        n = kk[n.index];
        
        pdf.name = paste("myplots/FIGURES/Targets_Kinases_", rownames(infer)[n],".pdf", sep = "")
        pdf(pdf.name, width=1.8, height=1.6)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        mm = match(rownames(infer)[n], colnames(mat.oc))
        nn = which((mat.oc[,mm])>0)
        targets = rownames(mat.oc)[nn]
        index = match(targets, phospho$seq.names)
        
        cbind(phospho$gene[index], phospho$phase[index], phospho$pval[index], rsum[nn])
        
        #index = index[which(phospho$qv[index]<0.05)]
        
        test = as.matrix(phospho[index, c(1:16)])
        length(index)
        test = t(apply(test, 1, standardization.nona))
        
        time = seq(0, 45, by=3)
        
        lims = range(test, na.rm=TRUE)
        
        ttest = unlist(strsplit(as.character(mtfs[n.index]), '_'))
        ttest = paste(ttest[-length(ttest)], sep='', collapse='_')
        
        matplot(c(0:15)*3, t(test), type='n', col="darkblue",lwd=1.2,cex=0.6, main=toupper(ttest), cex.main=0.8, ylim=lims, ylab=NA, xlab=NA, xlim=c(0,48), axes=FALSE)
        for(ii in 1:nrow(test))
        {
            points(time, test[ii, ], type='l', lwd=1.0, col='skyblue2')
        }
        
        ss = apply(test, 2, mean, na.rm=TRUE)
        points(time, ss, type='l', lwd=2.0, col='black')
        
        axis(1,at=seq(0, 48, by=12),cex.axis =1.0)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = signif(lims, d=1)
        if(gene=='Hdac3') {
            axis(2,at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),las=1,cex.axis = 1.0)
        }else{
            by = signif((lims[2]-lims[1])/5,d=0)
            axis(2,at = seq(lims[1], lims[2], by=by),las=1,cex.axis = 1.0)
        }
        box()
        abline(h=0,lty=2,lwd=1.5, col="darkgray")
        #if(gene=='Rorc')
        #legend(0,(lims[2]+0.5),c('phospho','nuclear'),pch=c(16,16),lty=c(1,1),col=c("darkgreen", "darkblue"), bty='n')
        
        dev.off()
        
    }
    
    
}

#################
####### FIGURE 6 (Main + SUpplementary)
#################
FIGURE_6 = TRUE
if(FIGURE_6)
{
    ############
    ##### clean the tables of TFs and coregulators
    ############
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
    
    #write.table(table.sx, file='Tables_DATA/Table_Supp_rhythmic_nuclear_proteins_wt_ko_test.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    
    ### rhythmic TFs detected in nuclear and phospho data
    cutoff.amp = 0.
    kk = which(!is.na(table.sx$TFs)==TRUE & table.sx$amp>=cutoff.amp)
    tfs = table.sx[kk,]
    kk = which(tfs$Gene.names=='Trp53;Tp53')
    if(length(kk)==1) tfs[kk, 2] = 'Tp53'
    
    Add.phospho = TRUE
    if(Add.phospho)
    {
        load(file='Rdata/Phos_Nuclear_Mapping_all_v2.Rdata')
        #load(file='Rdata/Phos_Nuclear_Mapping_all_v2.Rdata')
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
        #save(tfs.phospho, phospho, phospho.names, phospho.nuclear, file='Rdata/Phospho_Nuclear_Mapping_TFs_FDR_0.05.Rdata')
        load(file='Rdata/Phospho_Nuclear_Mapping_TFs_FDR_0.05.Rdata')
        
        tfs.names = unique(tfs.phospho$gene)
        mm = match(tfs.names, tfs$Gene.names)
        tfs.names[which(is.na(mm))]
        tfs.add = tfs.names[which(is.na(mm))]
        Add = matrix(NA, nrow=length(tfs.add), ncol=ncol(tfs))
        colnames(Add) = colnames(tfs)
        Add[,2] = tfs.add
        Add[, 10] = 1
        xx = data.frame(rbind(tfs, Add), stringsAsFactors=FALSE) ### tfs detected in nuclear and phospho proteins
        
        yy = matrix(NA, nrow=nrow(xx), ncol=3)
        colnames(yy) = c('amp.phospho', 'phase.phospho', 'qv.phospho')
        for(n in 1:nrow(yy))
        {
            gene = xx$Gene.names[n]
            ii = which(tfs.phospho$gene==gene)
            if(length(ii)==1)
            {
                yy[n,1] = tfs.phospho$amp[ii];
                yy[n,2] = tfs.phospho$phase[ii];
                yy[n,3] = tfs.phospho$qv[ii];
            }
            if(length(ii)>1)
            {
                print(gene);
                if(!is.na(xx$phase[n])==TRUE)
                {
                    diff =  tfs.phospho$phase[ii]-as.numeric(xx$phase[n]);
                    kk = which(diff>12)
                    diff[kk] = diff[kk]-12
                    kk = which(diff<=-12)
                    diff[kk] = diff[kk]+12
                    diff = abs(diff)
                    ii = ii[which(diff==min(diff))]
                    ii = ii[1]
                }else{
                    ii = ii[which(tfs.phospho$qv[ii]==min(tfs.phospho$qv[ii]))]
                    ii = ii[1]
                }
                
                yy[n,1] = tfs.phospho$amp[ii];
                yy[n,2] = tfs.phospho$phase[ii];
                yy[n,3] = tfs.phospho$qv[ii];
            }
        }
        
        xx = data.frame(xx, yy, stringsAsFactors=FALSE)
        
        kk = which(!is.na(xx$phase) & !is.na(xx$phase.phospho))
        cbind(xx$Gene.names[kk], xx$phase[kk], xx$phase.phospho[kk])
        
        tfs = xx
        
        length(which(!is.na(xx$phase) & is.na(xx$phase.phospho)))
        length(which(!is.na(xx$phase) & !is.na(xx$phase.phospho)))
        length(which(is.na(xx$phase) & !is.na(xx$phase.phospho)))
        length(which(is.na(xx$phase) & is.na(xx$phase.phospho)))
        
        kk = which(is.na(tfs$phase) & !is.na(tfs$phase.phospho))
        tfs$phase[kk] = as.numeric(tfs$phase.phospho[kk])
        tfs$amp[kk] = as.numeric(tfs$amp.phospho[kk])
        
        tfs$phase = as.numeric(tfs$phase)
        tfs$amp = as.numeric(tfs$amp)
        
    }
    
    
    #####rhythmic coregulators
    kk = which(!is.na(table.sx$Transcription.Cofactors)==TRUE & table.sx$amp>=cutoff.amp)
    cofactors = table.sx[kk,]
    kk = which(cofactors$Gene.names=='Med20;Gm20517')
    if(length(kk)==1) cofactors[kk,2] = 'Med20'
    kk = which(cofactors$Gene.names=='Kiaa1737/CIPC')
    if(length(kk)==1) cofactors[kk,2] = 'Cipc'
    
    kk = which(cofactors$Gene.names=='1600027N09Rik;Mrgbp')
    if(length(kk)==1) cofactors[kk,2] = 'Mrgbp'
    
    kk = which(cofactors$Gene.names=='Arid5b')
    if(length(kk)>0) cofactors = cofactors[-kk,]
    
    #### Add motif information
    Add.motif.infos = TRUE
    if(Add.motif.infos)
    {
        #write.table(res, file='Motifs_inferred_tfs_mapping_curate_alpha_0.1.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
        curating.table.motifs = FALSE
        if(curating.table.motifs)
        {
            xx = read.table(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/Motifs_inferred_tfs_mapping_curate_alpha_0.1.txt', sep='\t', header=TRUE, as.is = c(1,2,5))
            length(unique(xx[,1]))
            
            mm = match(tfs$Gene.names, xx$TFs.mapping)
            missing = tfs$Gene.names[which(is.na(mm))]
            print(missing)
            
            xx = xx[, c(1:5)]
            write.table(xx, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/Motifs_inferred_tfs_mapping_curate_alpha_0.1_manual.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
        }
        #mapping.all = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Manually_curation_TFs_pools.txt', header=TRUE, sep='\t')
        xx = read.table(file='Transcription_network/motif_analysis_encode/TFs_motifs_acvity_mapping_dhs_encode_alpha_0.05_matrix_occurrence_5kb2TSS_manual_v1.txt', sep='\t', header=TRUE, as.is = c(1,2,5))
        length(unique(xx[,1]))
        #yy = xx[order(-xx[,4]), ]
        
        mm = match(tfs$Gene.names, xx$TFs.mapping)
        missing = tfs$Gene.names[which(is.na(mm))]
        length(missing)
        length(which(!is.na(mm)==TRUE))
        #print(missing[order(missing)])
        
        mm = match(xx$TFs.mapping, tfs$Gene.names)
        
        xx$phases.tfs.mapping = tfs$phase[mm]
        xx$pval.tfs.mapping = tfs$pval[mm]
        xx$qv.tfs.mapping = tfs$qv[mm]
        motifs = xx
        
        #### rename the motif names and also modify the phase of motifs with the same name
        ms = unique(motifs[,1])
        for(n in 1:length(ms))
        {
            mm = ms[n]
            #cat(mm, '\n')
            kk = which(motifs[,1]==mm)
            test = motifs[kk, 2]
            test = test[which(!is.na(test)==TRUE)]
            test = unique(test)
            if(length(test)>0) motifs[kk, 2] = test;
        }
        
        kk = which(!is.na(motifs$phases.tfs.mapping))
        motifs = rbind(motifs[kk,], motifs[-kk,])
        motifs[c(1:length(kk)),]
        
        kk = which(is.na(motifs$Motif.names.modif))
        motifs$Motif.names.modif[kk] = motifs$Motifs.inferred[kk]
        
        mm = match(tfs$Gene.names, motifs$TFs.mapping)
        ff = motifs[mm, 2] 
        
        #### Check again the mapping motifs for rhythmic TFs
        ii = which(is.na(tfs$phase) & !is.na(tfs$phase.phospho))
        tfs$Motifs = ff

    }
    
    ##########
    #### Plot motif activity
    ##########
    PLOT.MOTIF.Activity = TRUE
    if(PLOT.MOTIF.Activity)
    {
        xx = motifs
        motif.names = unique(motifs[,2])
        
        kk = match(motif.names, xx[,2])
        #motif.names = mmotif.names
        aa = motifs[kk, ]
        phase.m = aa[,3]
        amp.m = aa[,4]
        
        motif.amp = 4.0;
        
        o1 = order(phase.m)
        phase.m = phase.m[o1]
        amp.m = amp.m[o1]
        xx = phase.m
        yy = phase.m
        for(n in 1:(length(xx)-1))
        {
            diff = (xx[n+1]-xx[n]);
            if(diff>=24) diff = diff-24;
            if(diff<0.15) xx[n+1] = xx[n]+0.15;
            if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
        }
        phase.m = xx
        
        aa = aa[o1, ]
        
        motif.names = aa[,2]
        ### modify motif names
        Modify.motif.names = TRUE
        if(Modify.motif.names)
        {
            for(n in 1:length(motif.names))
            {
                test = unlist(strsplit(as.character(motif.names[n]), '_'))
                n.end = length(test)
                if(test[n.end]=='p2'|test[n.end]=='p3'|test[n.end]=='f1')
                {
                    test = paste(paste(test[-n.end], sep='', collapse='_'), sep='')
                    print(test)
                    motif.names[n] = test
                }
                
                test = unlist(strsplit(as.character(motif.names[n]), '[.]'))
                n.end = length(test)
                if(test[n.end]=='m')
                {
                    test = paste(paste(test[-n.end], sep='', collapse='.'), sep='')
                    print(test)
                    motif.names[n] = test
                }
            }

        }

        
        amp = motif.amp
        motif.a = amp*cos(2*pi/24*phase.m)
        motif.b = amp*sin(2*pi/24*phase.m)
        CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
        motif.aa = Re(CC)
        motif.bb = Im(CC)
        
        amp = motif.amp + 0.1
        motif.a = amp*cos(2*pi/24*phase.m)
        motif.b = amp*sin(2*pi/24*phase.m)
        CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
        motif.txt1.aa = Re(CC)
        motif.txt1.bb = Im(CC)
        
        amp = motif.amp-0.15
        motif.a = amp*cos(2*pi/24*phase.m)
        motif.b = amp*sin(2*pi/24*phase.m)
        CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
        motif.txt2.aa = Re(CC)
        motif.txt2.bb = Im(CC)
        
        rmm=max(abs(c(motif.aa,motif.bb)))+3.5
        rr=c(-rmm,rmm)
        xlim = rr
        ylim = rr
        
        pdf('myplots/FIGURES/Motif_activities_inferred_intron_elastic_net_alpha_0.05_curated_encode.pdf', width=5, height=5)
        par(cex = 0.7, las = 1, mgp = c(0.1,0.1,0), mar = c(0.1,0.1,0.1,0.1)+0.1, tcl = -0.3)
        
        plot(motif.aa, motif.bb, main=NA, type='n', xlim=xlim, ylim=ylim, axes=F, xlab='', ylab='', lwd=2, pch=23, cex=2.0)
        
        phi=seq(0,2*pi,len=1000)
        motif.cycle = motif.amp-0.15
        lines(motif.cycle*cos(phi), motif.cycle*sin(phi), col='darkgray', lwd=2.0)
        lines(motif.cycle*cos(phi)/5*4, motif.cycle*sin(phi)/5*4, col='darkgray', lwd=1., lty=2)
        lines(motif.cycle*cos(phi)/5*3, motif.cycle*sin(phi)/5*3, col='darkgray', lwd=1., lty=2)
        lines(motif.cycle*cos(phi)/5*2, motif.cycle*sin(phi)/5*2, col='darkgray', lwd=1., lty=2)
        lines(motif.cycle*cos(phi)/5*1, motif.cycle*sin(phi)/5*1, col='darkgray', lwd=1., lty=2)
        #text(x=, y=, labels=c('10', '8', '6', '4*10^-2'))
        dist = 0.2
        arrows((motif.amp-dist)*(-1), 0, (motif.amp-dist), 0,  col='darkgray', lty=2, lwd=1., length=0.0);
        arrows(0, (motif.amp-dist)*(-1), 0, (motif.amp-dist), col='darkgray', lty=2, lwd=1., length=0.0);
        #rainbow = rainbow(length(motif.names),s = 0.85, v = 0.85)
        for(n in 1:length(motif.names))
        {
            points(motif.aa[n], motif.bb[n], pch=21, cex=0.8, col='blue',bg='white')
            
            if(phase.m[n]<=6) srt = 90-phase.m[n]/24*360;
            if(phase.m[n]>6 & phase.m[n]<=12) srt = 90-phase.m[n]/24*360;
            if(phase.m[n]>12 & phase.m[n]<=18) srt = 270-phase.m[n]/24*360;
            if(phase.m[n]>=18) srt = 270-phase.m[n]/24*360;
            cex = 0.5
            col = 'black'
            if(phase.m[n]<=12)
            {
                adj = 0;
                text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=cex, col='darkblue', srt=srt, adj=adj)
                #adj = 1;
                #text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=cex, col=col, srt=srt, adj=adj)
            }else{
                adj = 1;
                text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=cex, col='darkblue', srt=srt, adj=adj)
                #text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=cex, col=col, srt=srt, adj=adj)
            }
            
            Add.amps = TRUE
            if(Add.amps)
            {
                ### add amplitudes of motifs (average amplitudes per motif)
                #amp = ((6+log10(amp.m[n]))/5)*motif.cycle
                amp = (amp.m[n]/0.1)*motif.cycle
                amp.a = amp*cos(2*pi/24*phase.m[n])
                amp.b = amp*sin(2*pi/24*phase.m[n])
                AA= (amp.a -1i*amp.b) * exp(1i*pi/2)
                amp.a = Re(AA)
                amp.b = Im(AA)
                lines(c(0, amp.a), c(0, amp.b), type='l', lwd=1.0, col='darkblue')
                
            }
        }
        
        dev.off()
        
    }
    
    ###### Plot TFs and cofactors
    #### configuration for the plot
    phase.tf = as.numeric(tfs$phase)
    o1 = order(phase.tf)
    phase.tf = phase.tf[o1]
    tfs = tfs[o1,]
    tfs.names = tfs$Gene.names
    ### modify the phase gaps in order to better represent TFs
    xx = phase.tf
    yy = phase.tf
    for(n in 1:(length(xx)-1))
    {
        diff = (xx[n+1]-xx[n]);
        #if(diff<0) diff = diff+24;
        if(diff>=24) diff = diff-24;
        if(diff<0.18) xx[n+1] = xx[n]+0.18;
        if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
    }
    phase.tf = xx
    
    tf.amp = 7.0
    tf.aa = tf.amp*cos(2*pi/24*phase.tf)
    tf.bb = tf.amp*sin(2*pi/24*phase.tf)
    CC= (tf.aa -1i*tf.bb) * exp(1i*pi/2)
    tf.a = Re(CC)
    tf.b = Im(CC)
    tfs.coord = cbind(tf.a, tf.b)
    rownames(tfs.coord) = tfs.names
    
    tf.cycle  = tf.amp-0.2
    tf.text = tf.amp+0.25
    tf.ko = tf.amp - 0.5
    
    tf.aa.text = tf.text*cos(2*pi/24*phase.tf)
    tf.bb.text = tf.text*sin(2*pi/24*phase.tf)
    CC= (tf.aa.text -1i*tf.bb.text) * exp(1i*pi/2)
    tf.aa.text = Re(CC)
    tf.bb.text = Im(CC)
    
    tf.aa.ko = tf.ko*cos(2*pi/24*phase.tf)
    tf.bb.ko = tf.ko*sin(2*pi/24*phase.tf)
    CC= (tf.aa.ko -1i*tf.bb.ko) * exp(1i*pi/2)
    tf.aa.ko = Re(CC)
    tf.bb.ko = Im(CC)
    
    ### rhythmic cofactors
    phase.cof = cofactors$phase
    
    o1 = order(phase.cof)
    phase.cof = phase.cof[o1]
    cofactors = cofactors[o1,]
    cof.names = cofactors$Gene.names
    tfs.cofs = read.table('Tables_DATA/TFs_Coregulators_Mapping.txt', sep='\t', header=TRUE)
    nn = match(cof.names, tfs.cofs[,1])
    tfs.cofs = tfs.cofs[nn,]
    ### modify the phase gaps in order to better present cofactors
    xx = phase.cof
    yy = phase.cof
    for(n in 1:(length(xx)-1))
    {
        diff = (xx[n+1]-xx[n]);
        #if(diff<0) diff = diff+24;
        if(diff>=24) diff = diff-24;
        if(diff<0.12) xx[n+1] = xx[n]+0.12;
        if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
    }
    phase.cof = xx
    
    cof.amp = 10.0
    cof.a = cof.amp*cos(2*pi/24*phase.cof)
    cof.b = cof.amp*sin(2*pi/24*phase.cof)
    CC= (cof.a -1i*cof.b) * exp(1i*pi/2)
    cof.a = Re(CC)
    cof.b = Im(CC)
    
    cof.text = cof.amp+0.25
    cof.cycle  = cof.amp-0.2
    cof.ko = cof.amp - 0.5
    
    cof.aa.text = cof.text*cos(2*pi/24*phase.cof)
    cof.bb.text = cof.text*sin(2*pi/24*phase.cof)
    CC= (cof.aa.text -1i*cof.bb.text) * exp(1i*pi/2)
    cof.aa.text = Re(CC)
    cof.bb.text = Im(CC)
    
    cof.aa.ko = cof.ko*cos(2*pi/24*phase.cof)
    cof.bb.ko = cof.ko*sin(2*pi/24*phase.cof)
    CC= (cof.aa.ko -1i*cof.bb.ko) * exp(1i*pi/2)
    cof.aa.ko = Re(CC)
    cof.bb.ko = Im(CC)
    
    rmm=max(abs(c(tf.amp, tf.text, tf.cycle, cof.amp, cof.text, cof.cycle )))+1.5
    rr=c(-rmm,rmm)
    
    TFs.only = FALSE
    if(!TFs.only)
    {
        pdf.name = paste("myplots/FIGURES/TFs_Cofactors_v6.pdf")
        pdf(pdf.name, width=6.8, height=6.6)
        par(cex = 0.7, las = 1, mgp = c(0.1,0.1,0.1), mar = c(1.0, 1.0, 3.0, 3.0)+0.1, tcl = -0.3)
    }else{
        pdf.name = paste("myplots/FIGURES/TFs_v3.pdf")
        pdf(pdf.name, width=6.0, height=6.0)
        par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    }
    
    Add.phospho = TRUE
    Add.motif.infos = TRUE
    if(Add.phospho & Add.motif.infos)
    {
        plot(tf.a, tf.b, main=NA, type='n', xlim=rr, ylim=rr, axes=FALSE, xlab='', ylab='', cex=2.0, lwd=1.5, pch=22, col='black', bg='magenta')
        
        for(n in 1:nrow(tfs))
        {
            col.tf = 'black';
            bg.tf = 'magenta'
            pch.tf = 22
            cex.tf = 1.2
            #if(is.na(tfs$Motifs[n])==TRUE)
            #{
            #   pch.tf = 24
            #    cex.tf = 1.7
            #}
            
            #if(!is.na(tfs$qv[n])==TRUE & !is.na(tfs$qv.phospho[n])==TRUE) bg.tf = 'magenta4'
            if(!is.na(tfs$qv[n])==TRUE & !is.na(tfs$qv.phospho[n])==TRUE) bg.tf = 'slateblue4'
            #if(!is.na(tfs$qv[n])==FALSE & !is.na(tfs$qv.phospho[n])==TRUE) bg.tf = 'mistyrose2'
            if(!is.na(tfs$qv[n])==FALSE & !is.na(tfs$qv.phospho[n])==TRUE) bg.tf = 'yellow'
            #if(!is.na(tfs$relamp[n])==TRUE & !is.na(tfs$phase.phospho[n])==TRUE) col.tf = 'magenta'
            #if(is.na(tfs$relamp[n])==TRUE & !is.na(tfs$phase.phospho[n])==TRUE) col.tf = 'gray60'
            
            points(tf.a[n], tf.b[n], type='p', cex=cex.tf, pch=pch.tf, col=col.tf, bg=bg.tf)
            
            if(!is.na(tfs$prob.rhythmic.same.parameters.wt.ko[n])==TRUE)
            {
                if(as.numeric(tfs$prob.rhythmic.same.parameters.wt.ko[n])>0.5)
                {
                    points(tf.aa.ko[n], tf.bb.ko[n], pch=8, cex=1., col='black')
                }
            }
            ### add amplitudes of TFs
            amp = tfs$amp[n]/max(tfs$amp)*tf.cycle
            amp.a = amp*cos(2*pi/24*phase.tf[n])
            amp.b = amp*sin(2*pi/24*phase.tf[n])
            AA= (amp.a -1i*amp.b) * exp(1i*pi/2)
            amp.a = Re(AA)
            amp.b = Im(AA)
            lines(c(0, amp.a), c(0, amp.b), type='l', lwd=1.5, col=bg.tf)
            
        }
        
    }else{
        plot(tf.a, tf.b, main=NA, type='p', xlim=rr, ylim=rr, axes=FALSE, xlab='', ylab='', cex=2.0, lwd=1.5, pch=22, col='black',bg='magenta')
    }
    
    dist = 0.2
    arrows((tf.cycle)*(-1), 0, (tf.cycle), 0,  col='darkgray', lty=2, lwd=1., length=0.0);
    arrows(0, (tf.cycle)*(-1), 0, (tf.cycle), col='darkgray', lty=2, lwd=1, length=0.0);
    
    lwd = 1.0
    phi=seq(0,2*pi,len=1000)
    lines(tf.cycle*cos(phi), tf.cycle*sin(phi), col='darkgray', lwd=2)
    lines(tf.cycle*4/5*cos(phi), 4/5*tf.cycle*sin(phi), col='darkgray', lwd=lwd,lty=2)
    lines(tf.cycle*3/5*cos(phi), 3/5*tf.cycle*sin(phi), col='darkgray', lwd=lwd,lty=2)
    lines(tf.cycle*2/5*cos(phi), 2/5*tf.cycle*sin(phi), col='darkgray', lwd=lwd,lty=2)
    lines(tf.cycle*1/5*cos(phi), 1/5*tf.cycle*sin(phi), col='darkgray', lwd=lwd,lty=2)
    
    col='darkblue'
    
    for(n in 1:nrow(tfs))
    {
        if(phase.tf[n]<=6) srt = 90-phase.tf[n]/24*360;
        if(phase.tf[n]>6 & phase.tf[n]<=12) srt = 90-phase.tf[n]/24*360;
        if(phase.tf[n]>12 & phase.tf[n]<=18) srt = 270-phase.tf[n]/24*360;
        if(phase.tf[n]>=18) srt = 270-phase.tf[n]/24*360;
        
        if(phase.tf[n]<=12.0) {
            adj = 0; pos=2;
        }else{
            adj = 1; pos=4;
        }
        #text(tf.aa.text[n], tf.bb.text[n], tfs.names[n], cex=1.3,col=col, srt=srt, adj=adj, pos=pos)
        text(tf.aa.text[n], tf.bb.text[n], toupper(tfs.names[n]), cex=0.75,col=col, srt=srt, adj=adj)
        
    }
    
    if(TFs.only)
    {
        legend('topright', legend = c('TF.pred.MA','TF.non.pred.MA') , cex=1.2, pch=c(22, 22), col='darkblue', pt.cex=2.0, pt.lwd=1.5, pt.bg= c('magenta1', 'plum'), border = NA, bty = 'n')
        dev.off()
        
    }else{
        cex.cof = 1.2
        lines(cof.cycle*cos(phi), cof.cycle*sin(phi), col='darkgray', lwd=2)
        for(n in 1:nrow(cofactors))
        {
            col.bg = c('white','darkgreen','lawngreen','midnightblue','lightblue1', 'orange4','orange')
            col.index = 1
            if(!is.na(cofactors$Enzyme.Activity[n])==TRUE)
            {
                if(cofactors$Enzyme.Activity[n]=='A') col.index = 2
                if(cofactors$Enzyme.Activity[n]=='DA') col.index = 3
                if(cofactors$Enzyme.Activity[n]=='M') col.index = 4
                if(cofactors$Enzyme.Activity[n]=='DM') col.index = 5
                if(cofactors$Enzyme.Activity[n]=='U') col.index = 6
                if(cofactors$Enzyme.Activity[n]=='DU') col.index = 7
                #if(cofactors$Enzyme.Activity[n]=='K') col.index = 8
                #if(cofactors$Enzyme.Activity[n]=='P') col.index = 9
            }
            points(cof.a[n], cof.b[n], main=NA, type='p', cex=cex.cof, pch=21, col='black', bg=col.bg[col.index])
            
            if(phase.cof[n]<=6) srt = 90-phase.cof[n]/24*360;
            if(phase.cof[n]>6 & phase.cof[n]<=12) srt = 90-phase.cof[n]/24*360;
            if(phase.cof[n]>12 & phase.cof[n]<=18) srt = 270-phase.cof[n]/24*360;
            if(phase.cof[n]>=18) srt = 270-phase.cof[n]/24*360;
            
            if(phase.cof[n]<=12.0) {
                adj = 0;
            }else{
                adj = 1;
            }
            #text(cof.aa.text[n], cof.bb.text[n], cof.names[n], cex=1.2,col=col, srt=srt, adj=0, pos=pos)
            text(cof.aa.text[n], cof.bb.text[n], toupper(cof.names[n]), cex=0.75, col=col, srt=srt, adj=adj)
            
            if(!is.na(cofactors$prob.rhythmic.same.parameters.wt.ko[n])==TRUE)
            {
                if(as.numeric(cofactors$prob.rhythmic.same.parameters.wt.ko[n])>0.5)
                {
                    points(cof.aa.ko[n], cof.bb.ko[n], pch=8, cex=1., col='black')
                }
            }
            #if(!is.na(cofactors$pval.cryko[n])==TRUE & !is.na(cofactors$cor.cryko[n]))
            #{if(as.numeric(cofactors$pval.cryko[n])>0.05 & as.numeric(cofactors$cor.cryko[n])>=0.5)
            #}
            
            
        }
        
        #legend(x=6.5,y=12.1, legend = c(' ',' ', ' ') , cex=1.2, pch=c(22, 22, 22), col='darkblue', pt.cex=2.0, pt.lwd=1.5, pt.bg= c('magenta', 'slateblue4', 'yellow'), border = NA, bty = 'n')
        cex = 1.0
        pt.lwd=1.5
        pt.cex = 1.5
        legend(x=6.,y=12.8, legend = c('TF', 'nuclear','both', 'phospho') , pch=c(22, 22, 22, 22), col='darkblue', cex=cex, pt.cex=pt.cex, pt.lwd=pt.lwd, pt.bg= c('white', 'magenta', 'slateblue4', 'yellow'), border = NA, bty = 'n')
        legend(x=9,y=12.8, legend = c('coregulator', 'acetyl', 'deacetyl', 'methyl', 'demethyl', 'ubiq', 'deubiq'), cex=cex, pt.cex=pt.cex, pt.lwd=pt.lwd, pch=c(rep(21,7)), col='darkblue',  pt.bg= c(col.bg), border = NA, bty = 'n')
        legend(x=9, y=7.8, legend='no diff in WT and KO', cex=cex, pch=8, col='black', border = NA, bty = 'n')
        #text(x=6.8, y=12.1, 'M', col='black', cex=1.2)
        #text(x=7.4, y=12.1, 'NM', col='black', cex=1.2)
        dev.off()
        
    }
    
    
    ############
    ##### Compare TFs phases and motif activities
    ###########
    
    ### find paires of TFs and associated motifs and annotation of activator and repressor
    by.Function = TRUE
    if(by.Function)
    {
        ### TFs with corresponding motifs
        kk = match(tfs$Gene.names, motifs$TFs.mapping)
        tfs$Motifs = motifs[kk, 2]
        
        jj = match(motifs$TFs.mapping, tfs$Gene.names)
        jj = which(!is.na(jj)==TRUE)
        yy = motifs[jj,]
        
        rad = read.table(file='Tables_DATA/Activity_TFs_with_motifs.txt', sep='\t',header=TRUE)
        mm = match(yy$TFs.mapping, rad[,1])
        yy$activity = rad[mm, 2]
        
        col.a = 'green'
        col.r = 'red'
        col.d = 'orange'
        
        #### Phase delay between TFs and motif activity
        Plot.delay.Motifs.TFs = FALSE
        if(Plot.delay.Motifs.TFs)
        {
            pp.diff = yy$phases.tfs.mapping - yy$Phase.motifs.inferred
            kk = which(pp.diff<=(-6))
            pp.diff[kk] = pp.diff[kk]+24
            kk = which(pp.diff>18)
            pp.diff[kk] = pp.diff[kk]-24
            
            cex.axis = 1.0
            ylim = c(0, 15)
            pdfname = paste('myplots/FIGURES/Phase_Diff_TFs_motifs_all.pdf', sep='')
            pdf(pdfname, width=1.8, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            hist(pp.diff, breaks=seq(-6, 18, by=3), col='white', ylim=ylim, xlab=NA, ylab=NA,main=NA, axes=FALSE)
            box()
            axis(1,at=seq(-6,18,by=6),cex.axis =cex.axis)
            axis(2,at=seq(0,15,by=5),las=1,cex.axis = cex.axis)
            #abline(v=18, col='darkgray', lwd=1.5, lty=2)
            #abline(v=6, col='darkgray', lwd=2.0, lty=1)
            
            dev.off()
            
            
            pdfname = paste('myplots/FIGURES/Phase_Diff_TFs_motifs_activator.pdf', sep='')
            pdf(pdfname, width=1.8, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            #ylim = c(0, 10)
            kk = which(yy$activity=='Activator')
            hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=col.a, ylim=ylim, xlab=NA, ylab=NA,main=NA, axes=FALSE)
            #kk = which(yy$activity=='Repressor')
            #hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=rgb(1,0,0,0.6), xlab=NA, ylab=NA,main=NA, add=TRUE)
            box()
            axis(1,at=seq(-6,18,by=6),cex.axis =cex.axis)
            axis(2,at=seq(0,15,by=5),las=1,cex.axis = cex.axis)
            #abline(v=18, col='darkgray', lwd=1.5, lty=2)
            abline(v=6, col='darkgray', lwd=2.0, lty=1)
            
            dev.off()
            
            pdfname = paste('myplots/FIGURES/Phase_Diff_TFs_motifs_repressor.pdf', sep='')
            pdf(pdfname, width=1.8, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #kk = which(yy$activity=='Activator')
            #hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=rgb(0,1,0,0.6), xlab=NA, ylab=NA,main=NA, axes=FALSE)
            kk = which(yy$activity=='Repressor')
            hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=col.r, ylim=ylim, xlab=NA, ylab=NA,main=NA, axes=FALSE)
            box()
            axis(1,at=seq(-6,18,by=6),cex.axis =cex.axis)
            axis(2,at=seq(0,15,by=5),las=1,cex.axis = cex.axis)
            #abline(v=18, col='darkgray', lwd=1.5, lty=2)
            abline(v=6, col='darkgray', lwd=2.0, lty=1)
            
            dev.off()

            
            pdfname = paste('myplots/FIGURES/Phase_Diff_TFs_motifs_dual.pdf', sep='')
            pdf(pdfname, width=1.5, height=1.5)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            kk = which(yy$activity=='Dual')
            hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=col.d, xlab=NA, ylab=NA,main=NA, axes=FALSE)
            #kk = which(yy$activity=='Repressor')
            #hist(pp.diff[kk], breaks=seq(-6, 18, by=3), col=rgb(1,0,0,0.6), xlab=NA, ylab=NA,main=NA, add=TRUE)
            box()
            axis(1,at=seq(-6,18,by=6),cex.axis =0.8)
            axis(2,at=seq(0,10,by=5),las=1,cex.axis = 0.8)
            #abline(v=18, col='darkgray', lwd=1.5, lty=2)
            abline(v=6, col='darkgray', lwd=2.0, lty=1)
            
            dev.off()
            
        }
        
        ### regroup TFs with corresponding motifs
        kk = match(tfs$Gene.names, motifs$TFs.mapping)
        tfs$Motifs = motifs[kk, 2]
        jj = match(motifs$TFs.mapping, tfs$Gene.names)
        jj = which(!is.na(jj)==TRUE)
        yy = motifs[jj,]
        rad = read.table(file='Tables_DATA/Activity_TFs_with_motifs.txt', sep='\t',header=TRUE)
        mm = match(yy$TFs.mapping, rad[,1])
        yy$activity = rad[mm, 2]
        tfs.motifs.sel = yy
        
        #motifs[kk, 2] = motifs[kk, 1]
        #jj = which(!is.na(tfs$Motifs)==TRUE & tfs$Gene.names!='Esr1' & tfs$Gene.names!='Esrra' & tfs$Gene.names!='Nr2f6')
        #ttfs = tfs[jj,]
        #ttfs = ttfs[order(ttfs$DNA.Binding.Domains),]
        
        #kk = which(is.na(motifs[,2]))
        #motifs[kk, 2] = motifs[kk, 1]
        #familly = unique(table.sx$DNA.Binding.Domains)
        #kk = match(c('Esr1', 'Esrra'), )
        
        ### annotate 6 groups of TFs with associated motifs
        g1 = c('Arntl', 'Clock', 'Nr1d1', 'Nr1d2', 'Rora', 'Rorc') ### core clock
        #mm1 = c("E-Box", "E-Box", "RRE", "RRE", "RRE", "RRE")
        #rad1 = c('A', 'A', 'R', 'R', 'A', 'A')
        g2 = c('Dbp', 'Tef', 'Hlf', 'Nfil3', 'Bhlhe40', 'Klf3', 'Klf13') ### clock output regulators
        #mm2 = c("D-Box", "D-Box", "D-Box", "D-Box", "E-Box", "KLFE", "KLFE")
        #g2 = c('Dbp', 'Tef', 'Hlf', 'Nfil3', 'Bhlhe40') ### clock output regulators
        #mm2 = c("D-Box", "D-Box", "D-Box", "D-Box", "E-Box")
        #rad2 = c()
        g3 = c('Nr3c1', 'Nr3c2', 'Hnf4a', 'Nr1h4', 'Rxra', 'Esr1', 'Esrra', 'Ppara', 'Ppard') ### nuclear receptor
        #mm3 = c("IR3","IR3", "Hnf4s", "Nr1h4","RXRE")
        
        #g3 = c('Nr3c1', 'Nr3c2') ### nuclear receptor
        #mm3 = c("IR3","IR3")
        
        g4 = c('Atf1', 'Hnf1a', 'Hnf1b', 'Tfeb', 'Hsf1', 'Tcf12')
        #mm5 = c("CRE", "Hnf1s","Hnf1s", "Tfeb","HSE","E-Box")

        g5 = c('Foxa3', 'Foxk1', 'Foxp1', 'Foxp2', 'Foxa1')
        #mm4 = c("FOX_D1_D2_p2", "FOXK", "FOX_D1_D2_p2","FOX_D1_D2_p2","FOX_D1_D2_p2")

        #g4 = c('Atf1', 'Tfeb', 'Hsf1', 'Tcf12')
        #mm5 = c("ATF6_p2","TFEB.m","HSF1_2_p2","E-Box")
        
        g6 = c('Erf', 'Elf1', 'Elf2', 'Etv3', 'Fli1')
        #mm6 = c("ETS", "ETS", "ETS", "ETS", "ETS")
        #setdiff(c(g1, g2, g3, g4, g5, g6), yy$TFs.mapping)
        
        #setdiff(ttfs$Gene.names, c(g1, g2, g3, g4, g5, g6))
        
        g7 = setdiff(yy$TFs.mapping, c(g1, g2, g3, g4, g5, g6))
        #mm7 = c("IKZF1_p2","LEF1_TCF7_TCF7L1_2_p2", "Arid5s", "MEF2_A_B_C_D_p2", "E-Box", "MAFB.m", "CEBPA_B_DDIT3_p2")
        gg = c(g1, g2, g3, g4, g5, g6, g7)
        
        #rad = (c)
        
        Add.unknown.TF.activity = FALSE
        if(Add.unknown.TF.activity)
        {
            #f1 = c('bHLH', 'bZIP', 'HSF')
            #f2 = c('nuclearreceptor', 'forkhead', 'ETS')
            #rad = read.table(file='Tables_DATA/Activity_TFs_with_motifs.txt', sep='\t',header=TRUE)
            #mm = match(ttfs$Gene.names, rad[,1])
            #ttfs = data.frame(ttfs, rad[mm], stringsAsFactors=FALSE)
            #ttfs$activity = rad[mm, 2]
            #mm = match(ttfs$Gene.names, motifs$TFs.mapping)
            #ttfs = data.frame(ttfs, motifs[mm, c(2:5)], stringsAsFactors=FALSE)
            
            ## add unknown annotation
            yy = tfs.motifs.sel
            keep = yy$activity
            keep = as.character(keep)
            
            for(n in 1:nrow(yy))
            {
                pp.diff = as.numeric(yy$phases.tfs.mapping[n]) - as.numeric(yy$Phase.motifs.inferred[n])
                if(pp.diff<(-6)) pp.diff = pp.diff+24
                if(pp.diff>18) pp.diff = pp.diff-24
                #print(pp.diff)
                
                if(!is.na(yy$activity[n])==TRUE)
                {
                    if(yy$activity[n]=='Activator' & pp.diff>6)
                    {
                        keep[n] = 'Unknown'
                        #cat(pp.diff, ttfs[n, 2], '\n');
                    }
                    if(yy$activity[n]=='Repressor' & pp.diff<=6)
                    {
                        
                        keep[n] = 'Unknown'
                        #cat(pp.diff, ttfs[n, 2], '\n');
                    }
                    
                }
            }

            yy$activity = keep
            tfs.motifs.sel = yy
            tfs.motifs.sel[which(tfs.motifs.sel$activity=='Unknown'), ]
            #ff4$Phase.motifs.inferred = mean(ff4$Phase.motifs.inferred)
        }
        
        for(nn in 1:7)
        {
            eval(parse(text = paste('gg = g', nn, sep='')))
            mm = match(gg, yy$TFs.mapping)
            ff = yy[mm,]
            #kk = match(ff$Gene.names, motifs$TFs.mapping)
            ff = data.frame(ff, stringsAsFactors=FALSE)
            #eval(parse(text = paste('ff$Motif.names.modif = mm', nn, sep='')))
            eval(parse(text = paste('ff', nn, ' = ff', sep='')))
        }

        ###### for version of the encode DHS
        Motifs.used = FALSE
        if(Motifs.used)
        {
            mmref = unique(c(mm1, mm2, mm3, mm4, mm5, mm6, mm7))
            keep = c()
            for(rr in mmref)
            {
               keep = c(keep, motifs[c(which(motifs[,2]==rr), grep(toupper(rr), motifs[,2])),1])

            }
            keep = unique(keep)
            motifs.used = keep
            motifs.used = c(motif.used, c('MEF2_A_B_C_D_p2', 'ARID5B_p2'))
            save(motifs.used, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/Rdata/motifs_used_ref.Rdata')
        }
        
    }
    
    Table.Sup = FALSE
    if(Table.Sup)
    {
        ### motifs
        xx = motifs
        xx = xx[, -1]
        colnames(xx) = c('motifs', 'phase.motifs', 'amp.motifs', 'TFs.associated', 'phase.TFs', 'pval.TFs', 'qv.TFs')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S6_motif_activity_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        ### TFs with corresponding motifs and cofactors
        xx = tfs
        kk = c(2:3, 4:9, 21:24)
        xx = xx[,kk]
        colnames(xx)[3:8] = paste(colnames(xx[3:8]), '.nuclear', sep='')
        #colnames(xx)[9:11] = c('chow.test.WT.KO', 'nb.timepoints.KO', 'cor.WT.KO')
        
        yy = rep(NA, nrow(xx))
        yy[which(!is.na(xx$relamp.nuclear) & !is.na(xx$phase.phospho))] = 'both'
        yy[which(!is.na(xx$relamp.nuclear) & is.na(xx$phase.phospho))] = 'nuclear'
        yy[which(is.na(xx$relamp.nuclear) & !is.na(xx$phase.phospho))] = 'phospho'
       
        jj = which(is.na(xx$relamp.nuclear) & !is.na(xx$phase.phospho))
        xx$phase.nuclear[jj] = NA
        xx$amp.nuclear[jj] = NA
        xx[,2] = yy
        colnames(xx)[2] = 'rhythmic.in.nuclear.or.phospho'
        xx = xx[, c(1:2, 4, 6, 10:12, 9)]
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S6_TFs_nuclear_phospho_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        #index = table.sx[,1]
        #yy = nuclear[index, ]
        #yy = yy[, c(1:16, 34:37)]
        #xx = data.frame(table.sx[, -c(1, 12, 13, 19, 24:29)], xx, stringsAsFactors=FALSE)
        #colnames(xx)[3:8] = paste(colnames(xx[3:8]), '.WT', sep='')
        #colnames(xx)[17:19] = c('chow.test.WT.KO', 'nb.timepoints.KO', 'cor.WT.KO')
        #colnames(xx)[1:3] = c('kinase.motif', 'phase.LM', 'ampl.LM')
        xx = cofactors
        kk = c(2, 3:9, 21)
        xx = xx[,kk]
        xx = xx[, c(1, 4, 6, 9)]
        #colnames(xx)[4:6] = c('chow.test.WT.KO', 'nb.timepoints.KO', 'cor.WT.KO')
        #colnames(xx)[1:3] = c('kinase.motif', 'phase.LM', 'ampl.LM')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S6_coregulators_nuclear_all.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        xx = tfs.motifs.sel
        kk = c(5:9, 2:4)
        xx = xx[,kk]
        colnames(xx) = c('TFs', 'phase.TFs', 'pval.TFs', 'qv.TFs', 'activity.TFs', 'motifs.associated', 'phase.motifs', 'amp.motifs')
        #colnames(xx)[1:3] = c('kinase.motif', 'phase.LM', 'ampl.LM')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S6_TFs_motif_activity.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
    
    ######## PLot TFs phases vs. motif activity
    col.a = 'green'
    col.r = 'red'
    col.d = 'orange'
    for(nn in c(1:7))
    {
        pdfname = paste('myplots/FIGURES/TFs_vs_Motif_activities_ff', nn, '.pdf', sep='')
        
        eval(parse(text = paste('ff = ff', nn, sep='')))
        pp.diff = ff$phases.tfs.mapping - ff$Phase.motifs.inferred
        pp.diff[which(pp.diff<(-6))] = pp.diff[which(pp.diff<(-6))]+24
        pp.diff[which(pp.diff>=18)] = pp.diff[which(pp.diff>=18)]-24
        ff$phase.diff = pp.diff
        
        o1 = order(ff$phases.tfs.mapping)
        ff = ff[o1, ]
        tf.a = ff$TFs.mapping
        phase.tf = ff$phases.tfs.mapping
        motif.a = ff$Motif.names.modif
        phase.m = ff$Phase.motifs.inferred
        
        xx = phase.m
        yy = motif.a
        yy = unique(yy)
        
        if(length(yy)>1)
        {
            ii = match(yy, motif.a)
            xx = xx[ii]
            xx.o = order(xx)
            xx = xx[xx.o]
            yy = yy[xx.o]
            for(n in 1:(length(xx)-1))
            {
                diff = (xx[n+1]-xx[n]);
                if(diff>=24) diff = diff-24;
                #if(diff<0.5 & diff>0.01) xx[n+1] = xx[n]+0.5;
                if(diff<1) xx[n+1] = xx[n]+1;
                if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
            }
            ii = match(motif.a, yy)
            xx = xx[ii]
        }
        phase.m = xx
        
        xx = phase.tf
        for(n in 1:(length(xx)-1))
        {
            diff = (xx[n+1]-xx[n]);
            #if(diff<0) diff = diff+24;
            if(diff>=24) diff = diff-24;
            if(diff<0.5) xx[n+1] = xx[n]+0.5;
            if(xx[n+1]>=24) xx[n+1] = xx[n+1]-24;
        }
        phase.tf = xx
        
        pp.diff = ff$phase.diff
        
        tf.amp = 4.5
        tf.aa = tf.amp*cos(2*pi/24*phase.tf)
        tf.bb = tf.amp*sin(2*pi/24*phase.tf)
        CC= (tf.aa -1i*tf.bb) * exp(1i*pi/2)
        tf.aa = Re(CC)
        tf.bb = Im(CC)
        tf.amp = tf.amp+0.4
        tf.cc = tf.amp*cos(2*pi/24*phase.tf)
        tf.dd = tf.amp*sin(2*pi/24*phase.tf)
        CC= (tf.cc -1i*tf.dd) * exp(1i*pi/2)
        tf.cc = Re(CC)
        tf.dd = Im(CC)
        
        xlim = c(-(tf.amp+2.5), (tf.amp+2.5))
        ylim = xlim
        
        motif.amp = 3.0;
        motif.aa = motif.amp*cos(2*pi/24*phase.m)
        motif.bb = motif.amp*sin(2*pi/24*phase.m)
        CC = (motif.aa -1i*motif.bb) * exp(1i*pi/2)
        motif.aa = Re(CC)
        motif.bb = Im(CC)
        motif.txt = motif.amp-0.3
        motif.cc = (motif.txt)*cos(2*pi/24*phase.m)
        motif.dd = (motif.txt)*sin(2*pi/24*phase.m)
        CC = (motif.cc -1i*motif.dd) * exp(1i*pi/2)
        motif.cc = Re(CC)
        motif.dd = Im(CC)
        motif.cycle = motif.amp+0.2
        motif.ee = (motif.cycle)*cos(2*pi/24*phase.m)
        motif.ff = (motif.cycle)*sin(2*pi/24*phase.m)
        CC = (motif.ee -1i*motif.ff) * exp(1i*pi/2)
        motif.ee = Re(CC)
        motif.ff = Im(CC)
        
        if(nn<5)
        {
            pdf(pdfname, width=2.0, height=2.0)
        }else{
            pdf(pdfname, width=1.6, height=1.6)
        }
        #if(nn<5)
        #if(nn==5) pdf(pdfname, width=3.0, height=3.0)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.2,1.2,1.2,1.2)+0.1, tcl = -0.3)
        
        plot(tf.aa, tf.bb, main=NA, type='n', xlim=xlim, ylim=ylim, axes=F, xlab='', ylab='', lwd=1.2, pch=23, col='black',bg='green')
        #abline(v=0,h=0, col='darkgray',lwd=2.0)
        phi=seq(0,2*pi,len=1000)
        #lines((tf.amp-0.2)*cos(phi), (tf.amp-0.2)*sin(phi), col='darkgray', lwd=1)
        lines((tf.amp-0.2)*cos(phi), (tf.amp-0.2)*sin(phi), col='darkgray', lwd=1.5)
        lines((motif.cycle)*cos(phi), (motif.cycle)*sin(phi), col='darkgray', lwd=1.5)
        
        arrows((motif.cycle)*(-1), 0, (motif.cycle), 0,  col='darkgray', lty=2, lwd=1.0, length=0.0);
        arrows(0, (motif.cycle)*(-1), 0, (motif.cycle), col='darkgray', lty=2, lwd=1.0, length=0.0);
        
        for(n in 1:length(tf.a))
        {
            col='darkblue'
            if(phase.tf[n]<=6) srt = 90-phase.tf[n]/24*360;
            if(phase.tf[n]>6 & phase.tf[n]<=12) srt = 90-phase.tf[n]/24*360;
            if(phase.tf[n]>12 & phase.tf[n]<=18) srt = 270-phase.tf[n]/24*360;
            if(phase.tf[n]>=18) srt = 270-phase.tf[n]/24*360;
            if(phase.tf[n]<12.0) {adj=0; pos=4;}else{adj=1; pos = 2;}
            
            if(nn<5)
            {
                text(tf.cc[n], tf.dd[n], toupper(tf.a[n]), cex=0.7, col=col, srt=srt, adj=adj)
            }else{
                text(tf.cc[n], tf.dd[n], toupper(tf.a[n]), cex=0.6, col=col, srt=srt, adj=adj)
            }
            
            if(!is.na(motif.a[n]))
            {
                if(ff$activity[n]=='Activator') bg.col = col.a
                if(ff$activity[n]=='Repressor') bg.col = col.r
                if(ff$activity[n]=='Dual') bg.col = col.d
                
                if(ff$activity[n]=='Unknown') bg.col = 'black'
                arrows(tf.aa[n], tf.bb[n], motif.aa[n], motif.bb[n], col=bg.col, lty=1, lwd=1., length=0.01);
                if(ff$activity[n]=='Unknown') bg.col = 'white'
                points(tf.aa[n], tf.bb[n], lwd=1., pch=23, col='darkblue',bg=bg.col,cex=1.0)
                
            }else{
                bg.col = 'white';
                points(tf.aa[n], tf.bb[n], lwd=1.0, pch=23, col='darkblue',bg=bg.col,cex=1.0)
            }
        }
        mmotif.a = unique(motif.a)
        mmotif.a = mmotif.a[which(!is.na(mmotif.a))]
        phase.mm = phase.m[match(mmotif.a, motif.a)]
        mmotif.aa = motif.aa[match(mmotif.a, motif.a)]
        mmotif.bb = motif.bb[match(mmotif.a, motif.a)]
        mmotif.cc = motif.cc[match(mmotif.a, motif.a)]
        mmotif.dd = motif.dd[match(mmotif.a, motif.a)]
        for(n in 1:length(mmotif.a))
        {
            points(mmotif.aa[n], mmotif.bb[n], pch=21, cex=1.0, col='blue',bg='white')
            if(phase.mm[n]<=6) srt = 90-phase.mm[n]/24*360;
            if(phase.mm[n]>6 & phase.mm[n]<=12) srt = 90-phase.mm[n]/24*360;
            if(phase.mm[n]>12 & phase.mm[n]<=18) srt = 270-phase.mm[n]/24*360;
            if(phase.mm[n]>=18) srt = 270-phase.mm[n]/24*360;
            
            #col = 'black'
            #if(phase.mm[n]<=12) {adj = 1; pos=1;}else{adj = 0;pos=1;}
            #text(mmotif.cc[n], mmotif.dd[n], mmotif.a[n], cex=0.6, col=col, srt=srt, adj=adj)
            
        }
        #if(nn==1)   legend('topright', legend = c('activator', 'repressor' ), cex=1.2, pch=c(22), col=col, pt.cex=1.5, pt.lwd=1.2, pt.bg= c('green', 'red'), border = NA, bty = 'n')
        #if(nn==1) legend('topright', legend = c('dual' ), cex=1.25, pch=c(22), col=col, pt.cex=1.5, pt.lwd=1.2, pt.bg= c('orange'), border = NA, bty = 'n')
        
        dev.off()
        
    }
    
    ##########
    #### Examples of histone modifiers and histone modifications
    ##########
    hmr = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/DATA_WB_2.csv', sep=',', header=TRUE)
    colnames(hmr) = toupper(colnames(hmr))
    source('functions_nuclear.R')
    
    require('plotrix')
    examples2 = c('Sirt7', 'Hdac3', 'Sin3a')
    
    match(examples2, nuclear.names[,3])
    match(examples2, phospho.names[,2])
    
    cols = 'plum'
    
    for(gene in examples2)
    {
        pdf.name = paste("myplots/FIGURES/Fig_6_histone_modifiers_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.6, height=1.4)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        n = which(nuclear$Gene.names==gene)
        y0 = as.numeric(nuclear[n, c(1:16)])
        y0 = 2^y0
        y0 = y0/mean(y0, na.rm=TRUE)
        
        kk = grep(toupper(gene), colnames(hmr))
        
        y1 = rbind(as.numeric(unlist(hmr[c(1:16), kk[1]])), as.numeric(unlist(hmr[c(1:16), kk[2]])))
        #y1 = rbind(y1[c(1:8)], y1[c(9:16)], y1[17:24], y1[25:32])
        for(ii in 1:2)
        {
            y1[ii, ] = y1[ii, ]/mean(y1[ii, ], na.rm=TRUE)
        }
        
        averg1 = c()
        err1 = c()
        for(ii in 1:16)
        {
            xx = y1[,ii]
            #xx = xx/mean(xx[which(!is.na(xx)==TRUE)])
            ll = length(which(is.na(xx)==TRUE))
            if(ll==0){averg1 = c(averg1, mean(xx)); err1=c(err1, sd(xx)/sqrt(2))}
            if(ll==1){averg1 = c(averg1, mean(xx[which(!is.na(xx)==TRUE)])); err1=c(err1, sd(xx[which(!is.na(xx)==TRUE)])/sqrt(length(xx[which(!is.na(xx)==TRUE)])))}
            if(ll==2){averg1 = c(averg1, NA); err1=c(err1, NA)}
            
        }
        
        time = c(0:15)*3
        lims = range(c(y0, averg1-err1, averg1+err1), na.rm=TRUE)
        
        #lims = range(c(y0[which(!is.na(y0)==TRUE)]))
        plot(time, y0, type='l', ylim=lims, xlim=c(0, 48), col='darkblue', lwd=1.2, pch=16, main=NA, cex.main=0.8, axes=FALSE, xlab=NA, ylab=NA)
        points(time, y0, type='p', cex=0.7, pch=16, col='darkblue')
        points(time, averg1, type='l', lwd=1.2, col=cols)
        points(time, averg1, type='p', cex=0.7, col=cols, pch=15)
        #arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col='darkblue', lwd=1.5)
        arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3, col=cols, lwd=1.2, cex=0.7)
        
        #axis(1,at=12*c(0:4),cex.axis =1.)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
        by = 0.1
        print(gene)
        print(lims)
        
        cex = 0.8
        #axis(2,at = seq(lims[1], lims[2],by=0.3),las=1,cex.axis = 1.2)
        if(gene=='Sirt7') axis(2,at = c(0.6, 1, 1.4),las=1,cex.axis = cex)
        if(gene=='Hdac3') axis(2,at = c(0.6, 1, 1.4),las=1,cex.axis = cex)
        if(gene=='Sin3a') axis(2,at = c(0.6, 1, 1.4),las=1,cex.axis = cex)
        box()
        
        abline(h=1,lty=2,lwd=1.5, col="darkgray")
        dev.off()
        
    }
    
    hmd = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/H3modif.csv', sep=';', header=TRUE)
    rownames(hmd) = hmd[,1]
    hmd = hmd[, -1]
    hmd = hmd[-1, ]
    
    time = c(0:31)*3
    matplot(time, t(hmd), type='b', col=c(1:7), pch=c(1:7))
    res = t(apply(hmd, 1, f24_R2_alt2, t=time))
    
    pdf(paste('myplots/FIGURES/Histone_activation_markers.pdf', sep=''), width=2.0, height=1.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,3,2,0.8)+0.1, tcl = -0.3)
    
    time = c(0:7)*3
    lims = c(0.4, 2.0)
    
    for (n in c(1:5))
    {
        xx = mean.err(hmd[n,])
        averg = xx[1,]
        err = xx[2, ]
        col = n
        pch = 1
        lty = 1
        if(n == 1)
        {
            plot(time, averg, type='l', ylim=lims, xlim=c(0,24), col=col, lwd=1.2, lty=lty, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
            arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=lty, pch=n,lwd=1.2)
        }else{
            points(time, averg, type='l', lwd=1.2, col=col, pch=pch, lty=lty)
            arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, pch=n,lty=lty,lwd=1.2)
        }
    }
    
    #axis(1,at=6*c(0:4),cex.axis =1.0)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2,at = c(0.5, 1, 1.5, 2),las=1,cex.axis = 1.0)
    
    box()
    abline(h=1,lty=2,lwd=1.5, col="darkgray")
    
    legend(8, (lims[2]+0.1), c('H3K27ac','H3K4me1', 'H3K4me3', 'H3K36me3', 'H3K9ac'), lty=c(1),col=c(1:5), cex=0.6, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()
    
    
    hmd = hmd[c(6:7), ]
    
    pdf(paste('myplots/FIGURES/Histone_repression_markers.pdf', sep=''), width=2.0, height=1.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    time = c(0:7)*3
    lims = c(0.4, 1.8)
    
    for (n in c(1:2))
    {
        xx = mean.err(hmd[n,])
        averg = xx[1,]
        err = xx[2, ]
        if(n==1) col=6
        if(n==2) col = 8
        pch = 1
        lty = 1
        if(n == 1)
        {
            plot(time, averg, type='l', ylim=lims, xlim=c(0,24), col=col, lwd=1.2, lty=lty, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
            arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=lty, pch=pch,lwd=1.2)
        }else{
            points(time, averg, type='l', lwd=1.2, col=col, pch=pch, lty=lty)
            arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, pch=pch,lty=lty,lwd=1.2)
        }
    }
    axis(1,at=6*c(0:4),cex.axis =1.)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2,at = c(0.5, 1, 1.5, 2),las=1,cex.axis = 1.0)
    
    box()
    abline(h=1,lty=2,lwd=2.0, col="darkgray")
    
    legend(8, (lims[2]+0.05), c('H3K9me2/3','H3K27me3'), lty=c(1), col=(c(6, 8)), cex=0.6, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()
    
    ##################
    ### HATs and HDACs
    ##################
    hdac = c(paste('Hdac', c(1:11), sep=''), paste('Sirt', c(1:7), sep=''))
    #hat = c('Kat2b', 'Kat2a', 'Taf5l', 'Elp3', 'Hat1', 'Hat2','Ep300', 'Crebbp', 'Ncoa1', 'Ncoa3', 'Kat6a', 'Kat6b','Clock', 'Taf1', 'Kat8', 'Kat7', 'Kat5')
    hat = c('Hat1','Kat2b', 'Kat2a', 'Crebbp', 'Ep300', 'Taf1', 'Kat5', 'Kat5a', 'Kat5b', 'Kat7', 'Kat8', 'Elp3', 'Gtf3c4', 'Ncoa1', 'Ncoa2', 'Ncoa3')
    kmt = c('Suv39h1', 'Suv39h2', 'Ehmt2', 'Ehmt1', 'Setdb1', 'Setdb2', 'Mll', 'Mll2', 'Mll3', 'Mll4', 'Mll5', 'Setd1a', 'Setd1b', 'Setd2', 'Nsd1', 'Smyd2', 'Dot1l', 'Setd8', 'Suv420h1', 'Suv420h2', 'E2h2', 'Setd7', 'Prdm2')
    kdm = c('Kdm1a', 'Kdm1b', 'Kdm2a', 'Kdm2b', 'Kdm3a', 'Kdm3b', 'Kdm4a', 'Kdm4b', 'Kdm4c', 'Kdm5a', 'Kdm5b', 'Kdm5c', 'Kdm5d', 'Kdm6a', 'Kdm6b', 'Jhdm1d', 'Jmjd1c')
    
    enzymes = c('hdac', 'hat', 'kmt', 'kdm')
    for(enz in enzymes)
    {
        index = c()
        #eval(parse(text = paste('index', '=c()', sep='')))
        eval(parse(text = paste('xx = ', enz, sep='')))
        for(x in xx)
        {
            mm = nuclear.names[which(nuclear.names[,3]==x), 1]
            #print(mm)
            if(length(mm)>1) index = c(index, mm[which(nuclear$pval[mm]==min(nuclear$pval[mm]))])
            if(length(mm)==0) index = c(index, NA)
            if(length(mm)==1) index = c(index, mm)
        }
        
        index = index[which(!is.na(index)==TRUE)]
        index = index[which(nuclear$pval[index]<0.01)]
        eval(parse(text = paste('index.', enz, '=index', sep='')))
        #xx[which(is.na(index))]
        #nuclear[index,-c(1:19)]
    }
    
    col.bg = c('white','darkgreen','lawngreen','midnightblue','lightblue1', 'orange4','orange')
    col.bg = col.bg[c(2:5)]
    
    pdf('myplots/FIGURES/HDACs_HATs_phases.pdf', width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    ylim = range(nuclear$amp[c(index.hat, index.hdac)], na.rm=TRUE)
    ylim = c(0, 1.4)
    plot(nuclear$phase[index.hat], nuclear$amp[index.hat], xlim=c(0,24), col='darkblue', bg=col.bg[1], cex=1., pch=21, ylim = ylim, xlab=NA, ylab=NA, axes=FALSE, log='')
    points(nuclear$phase[index.hdac], nuclear$amp[index.hdac], col='darkblue', bg=col.bg[2], cex=1., pch=21)
    #abline(h=-log10(0.05), lwd=1.0, col='darkgray')
    
    #pos = 1;
    offset.init = 0.3
    cex = 0.6
    for(n in 1:length(index.hat))
    {
        index = index.hat[n]
        offset = offset.init;
        pos = 1;
        if(nuclear$Gene.names[index]=='Crebbp') pos=1
        if(nuclear$Gene.names[index]=='Ep300') pos=4;
        if(nuclear$Gene.names[index]=='Kat5') pos=2;
        text(nuclear$phase[index], nuclear$amp[index], toupper(nuclear$Gene.names[index]), cex=cex, col='black', pos=pos, offset=offset)
    }
    for(n in 1:length(index.hdac))
    {
        index = index.hdac[n]
        pos = 3;
        offset = offset.init
        if(nuclear$Gene.names[index]=='Hdac1;Gm10093')
        {
            pos=4;
            text(nuclear$phase[index], nuclear$amp[index], toupper('Hdac1'), cex=cex, col='black', pos=pos, offset=offset)
        }else{
            if(nuclear$Gene.names[index]=='Hdac5') pos=2;
            if(nuclear$Gene.names[index]=='Hdac3') pos=3;
            if(nuclear$Gene.names[index]=='Hdac2') pos=2;
            
            if(nuclear$Gene.names[index]=='Sirt2') pos=4;
            if(nuclear$Gene.names[index]=='Sirt7') pos=2;
            text(nuclear$phase[index], nuclear$amp[index], toupper(nuclear$Gene.names[index]), cex=cex, col='black', pos=pos, offset=offset)
        }
    }
    
    axis(1,at=6*c(0:4),cex.axis =1.)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    #lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2,at = seq(0,1.4, by=0.4),las=1,cex.axis = 1.)
    box()
    
    #abline(h=-log10(0.1), lwd=2.0, col='darkgray')
    dev.off()
    
    pdf('myplots/FIGURES/KMTs_KDMs_phases.pdf', width=1.8, height=1.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    ylim = range(nuclear$amp[c(index.kmt, index.kdm)], na.rm=TRUE)
    ylim = c(0.2, 0.65)
    plot(nuclear$phase[index.kmt], nuclear$amp[index.kmt], xlim=c(0,24), bg=col.bg[3], col='darkblue', cex=1., pch=21,  ylim = ylim, xlab=NA, ylab=NA, axes=FALSE)
    points(nuclear$phase[index.kdm], nuclear$amp[index.kdm], bg=col.bg[4],col='darkblue', cex=1., pch=21)
    
    offset.init = 0.3;
    pos = 1;
    cex = 0.6
    for(n in 1:length(index.kmt))
    {
        index = index.kmt[n]
        print(kmt[n])
        pos = 1;
        offset = offset.init
        #if(nuclear$Gene.names[index]=='Dot1l') pos=3;
        if(nuclear$Gene.names[index]=='Setdb1') pos=2;
        if(nuclear$Gene.names[index]=='Ehmt2') pos=4;
        if(nuclear$Gene.names[index]=='Setd1a') pos=4;
        if(nuclear$Gene.names[index]=='Setd1b') pos=1;
        if(nuclear$Gene.names[index]=='Setd8') pos=3;
        #if(nuclear$Gene.names[index]=='Mll5') pos=4;
        text(nuclear$phase[index], nuclear$amp[index], toupper(nuclear$Gene.names[index]), cex=cex, col='black', pos=pos, offset=offset)
    }
    for(n in 1:length(index.kdm))
    {
        index = index.kdm[n]
        print(kdm[n])
        pos = 3;
        offset = offset.init
        if(nuclear$Gene.names[index]=='Kdm5a') pos=3;
        if(nuclear$Gene.names[index]=='Kdm3a') pos=1;
        if(nuclear$Gene.names[index]=='Kdm1b') pos=2;
        if(nuclear$Gene.names[index]=='Jmjd1c') pos=2;
        if(nuclear$Gene.names[index]=='Kdm2a') pos=3;
        if(nuclear$Gene.names[index]=='Kdm6a') pos=4;
        text(nuclear$phase[index], nuclear$amp[index], toupper(nuclear$Gene.names[index]), cex=cex, col='black', pos=pos, offset=offset)
    }
    abline(h=-log10(0.05), lwd=1.0, col='darkgray')
    
    axis(1,at=6*c(0:4),cex.axis =1.)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    #lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2,at = seq(0,1.2, by=0.2),las=1,cex.axis = 1.)
    box()
    #text(nuclear$phase[index.kmt], -log10(nuclear$pval[index.kmt]), kmt, cex=0.7, col='black', pos=1, offset=0.3)
    #text(nuclear$phase[index.kdm], -log10(nuclear$pval[index.kdm]), kdm, cex=0.7, col='black', pos=3, offset=0.3)
    #abline(h=-log10(0.05), lwd=1.0, col='darkgray')
    #abline(h=-log10(0.1), lwd=2.0, col='darkgray')
    dev.off()
    
    
    #########
    #### Individual examples of TFs and cofactors for WT and KO
    #########
    cutoff.correlation = 0.5
    cutoff.pval = 0.05
    kk = which(tfs$pval.cryko<cutoff.pval & tfs$nb.cryko==4 & tfs$cor.cryko<cutoff.correlation)
    jj = which(cofactors$pval.cryko<cutoff.pval & cofactors$nb.cryko==4 & cofactors$cor.cryko<cutoff.correlation)
    
    ### individual examples
    require('plotrix')
    #examples = c('Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Nampt', 'Parp1', 'Prkaa2', 'Prkab1;Prkab2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d;Csnk1e')
    #cutoff.correlation = 0.7
    #kk = which(tfs$pval.cryko>0.1 & tfs$nb.ko==4 & tfs$cor.wtko>cutoff.correlation)
    #jj = which(cofactors$pval.cryko>0.1 & cofactors$nb.ko==4 & cofactors$cor.wtko>cutoff.correlation)
    #kk = kk[order(-tfs$cor.wtko[kk])]
    #jj = jj[order(-tfs$cor.wtko[jj])]
    
    #examples = c(nuclear$Gene.names[as.numeric(tfs[kk, 1])])
    #examples = c(nuclear$Gene.names[as.numeric(cofactors[jj,1])])
    #examples = examples[which(!is.na(examples)==TRUE)]
    #examples = setdiff(examples, c('Arntl', 'Clock', 'Per1', 'Per2', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp',  'Tef', 'Hlf', 'Nfil3'))
    examples = c('Bhlhe40', 'Ppara', 'Tcf12', 'E4f1', 'Mef2d', 'Srebf1', 'Elf1', 'Foxp1', 'Klf3', 'Nr2c1', 'Ehmt2', 'Hdac3', 'Nono', 'Sirt7', 'Trim24', 'Crebbp', 'Crtc2', 'Mll3', 'Epc2', 'Cxxc5')
    mm = match(examples, nuclear$Gene.names)
    
    ## remove examples in the following folder
    cmd = "rm -r /Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/myplots/FIGURES/TFs_cofactors_examples/*"
    system(cmd)
    
    for(gene in examples)
    {
        pdf.name = paste("myplots/FIGURES/TFs_cofactors_examples/TFs_Coregulators_examples_WT_KO_not_centered_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.4, height=1.1)
        par(cex = 0.7, las = 1, mgp = c(0.3,0.3,0), mar = c(1.7,1.7,1.0,0.8)+0.1, tcl = -0.3)
        
        n = which(nuclear$Gene.names==gene)
        if(length(n)>1) n = n[which(n==min(n))]
        name = nuclear[n,20]
        y01 = as.numeric(nuclear[n,c(1:8)])
        y02 = as.numeric(nuclear[n,c(9:16)])
        y0 = rbind(y01, y02)
        
        averg = c()
        err = c()
        for(ii in 1:8)
        {
            xx = y0[,ii]
            kk = length(which(is.na(xx)==TRUE))
            if(kk==0){averg = c(averg, mean(xx)); err=c(err, sd(xx)/sqrt(2))}
            if(kk==1){averg = c(averg, xx[which(!is.na(xx)==TRUE)]); err=c(err, 0)}
            if(kk==2){averg = c(averg, NA); err=c(err, NA)}
            
        }
        kk = grep('Cry.KO', colnames(nuclear))
        y1 = as.numeric(nuclear[n, kk])
        time = c(0:7)*3
        lims = range(c(averg-err, averg+err, y1), na.rm=TRUE)
        
        plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=toupper(gene), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
        points(time, averg, col="darkblue", lwd=1.2, type='l')
        points(6*c(0:3), y1, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
        cex = 0.8
        axis(1,at=6*c(0:4),cex.axis =cex)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = signif(lims, d=1)
        if(gene=='Srebf1') axis(2,at = c(-2, -1.5, -1, -0.5, 0, 0.5),las=1,cex.axis = cex)
        if(gene=='Elf1') axis(2,at = c(-0.4, -0.2, 0, 0.2),las=1,cex.axis = cex)
        if(gene=='Foxp1') axis(2,at = c(0.2, 0.4, 0.6),las=1,cex.axis = cex)
        if(gene=='Klf3') axis(2,at = c(-0.2, 0, 0.2),las=1,cex.axis = cex)
        if(gene=='Nr2c1') axis(2,at = c(-0.6, -0.4, -0.2, 0, 0.2),las=1,cex.axis = cex)
        if(gene=='Crebbp') axis(2,at = c(-0.2, -0.1, 0, 0.1),las=1,cex.axis = cex)
        if(gene=='Crtc2') axis(2,at = c(-0.6, -0.3, 0, 0.3, 0.6),las=1,cex.axis = cex)
        if(gene=='Bhlhe40') axis(2,at = c(-2.0, -1., 0, 1.0, 2.0),las=1,cex.axis = cex)
        if(gene=='Ppara') axis(2,at = c(-1.2,-0.8, -0.4, -0.4, 0, 0.4),las=1,cex.axis = cex)
        if(gene=='E4f1') axis(2,at = c(-0.4, 0, 0.4, 0.8),las=1,cex.axis = cex)
        if(gene=='Nono') axis(2,at = c(0, 0.1, 0.2),las=1,cex.axis = cex)
        if(gene=='Hdac3') axis(2,at = c(0, 0.2, 0.3, 0.4, 0.5, 0.6),las=1,cex.axis = cex)
        if(gene=='Sirt7') axis(2,at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),las=1,cex.axis = cex)
        if(gene=='Nr2c1'|gene=='Mll3'|gene=='Epc2'|gene=='Cxxc5') axis(2,at = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),las=1,cex.axis = cex)
        if(gene=='Tcf12'|gene=='Mef2d') axis(2,at = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),las=1,cex.axis = cex)
        if(gene=='Ehmt2'|gene=='Trim24') axis(2,at = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0),las=1,cex.axis = cex)
        
        box()
        
        abline(h=0,lty=2,lwd=1.5, col="darkgray")
        
        dev.off()
    }
    
    #### plot some examples of TFs and coregulators
    examples = c('Ppara', 'Ppard', 'Hdac3', 'Sirt7', 'Hsf1', 'Srebf1', 'Tfeb', 'Zkscan3')
    jj = match(examples, nuclear$Gene.names)
    
    for(n in jj)
    {
        gene = nuclear$Gene.names[n]
        
        print(gene)
        pdf.name = paste("myplots/FIGURES/TFs_coregulators_examples_Nuclear_centered_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.8, height=1.6)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        y0 = as.numeric(nuclear[n, c(1:16)])
        y0 = y0-mean(y0, na.rm=TRUE)
       
        time = c(0:15)*3
        lims = range(c(y0), na.rm=TRUE)
        
        xlim = c(0, 48)
        col = 'darkblue'
        pch = 16
        lwd = 1.2
        cex = 1.0
        plot(time, y0, type='l', ylim=lims, xlim=xlim, col=col, lwd=lwd,  pch=pch, lty=1, axes=FALSE, xlab=NA, ylab=NA,  main=toupper(gene), cex.main=0.8)
        #arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
        points(time, y0, type='p', col=col, pch=pch, cex=0.7)
        col='indianred1'
        #points(c(0,6,12, 18), y2, type='b', col=col, lwd=lwd, lty=2, pch=pch)
        
        axis(1,at=c(0:4)*12,cex.axis =cex)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        #lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
        
        diff = lims[2]-lims[1]
        if(diff>=3. & diff<6.0) axis(2,at = seq(-5, 5, by=1),las=1,cex.axis = cex)
        if(diff>=1.0 & diff<3.) axis(2,at = seq(-5, 5, by=0.5),las=1,cex.axis = cex)
        if(diff<1) axis(2,at = seq(-5, 5, by=0.2),las=1,cex.axis = 1.)
        box()
        abline(h=0,lty=2,lwd=2.0, col="darkgray")
        
        dev.off()
        
    }
    
}

#################
####### FIGURE 7 (Main + SUpplementary)
#################
FIGURE_7 = TRUE
if(FIGURE_7)
{
    Cell.cycle.markers = FALSE
    if(Cell.cycle.markers)
    {
        #### Protein complexes
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
        load(file='Rdata/Analysis_res_Protein_Complexes_Detected_Rhythmic_Selected_manual_v3.Rdata')
        res.c = res.c[, -c(5,6)]
        kk = which(res.c$amp.pc<=2.1)
        res.c = res.c[kk,]
        aa = res.m
        xx = res.c
        o1 = order(-xx$amp.pc)
        xx = xx[o1,]
        yy = xx
        #indexx = c(1:nrow(res.c))
        pc.names = c('DNMT1-RB1-HDAC1-E2F1', 'APC/C')
        mm = match(pc.names, xx$names.pc)
        xx = yy[mm, ]
        pc.test = xx
        #pc.names = res.c[,1]
        #pc.cols = c(rep('magenta', 3), rep('darkblue', 5), rep('black', 2))
        
        for(n in 1:nrow(pc.test))
        {
            name = pc.test[n, 1]
            jj = which(res.m[,1]==name)
            
            subunits = unique(unlist(strsplit(as.character(res.m$subunits[jj]), ',')))
            index = c()
            for(ss in subunits)
            {
                index = c(index, nuclear.names[which(nuclear.names[,3]==ss) ,1])
            }
            index = unique(index[which(!is.na(index)==TRUE)])
            #kk = unique(c(match(res.m$Complex.id[jj], res.sel$Complex.id), match(res.m$Complex.name[jj], res.sel$Complex.name)))
            #ii = unique(kk[which(!is.na(kk)==TRUE)])
            #if(length(ii)!=length(kk)) print(n)
            #kk = ii
            pc.name = c(as.character(pc.test[n, 5]), as.character(pc.test[n, 2]))
            pc.name = unlist(strsplit(as.character(pc.name), ' '))
            pc.name = unlist(strsplit(as.character(pc.name), '/'))
            pc.name = unlist(strsplit(as.character(pc.name), '-'))
            pc.name = paste(pc.name, sep='', collapse='_')
            pc.name = paste(pc.name, '_', signif(pc.test$amp.pc[n],d=2), sep='', collapse='_')
            
            #genes = unlist(strsplit(as.character(res.sel$subunits.detected[kk]), ','))
            kk = which(nuclear$nb.timepoints[index]>=8)
            index = index[kk]
            #genes = genes[kk]
            
            if(length(index)>1 & pc.test$amp.pc[n]>0.05)
            {
                pdf.name = paste("myplots/FIGURES/cell_markers_KO/", pc.name, "_Rhythmic_Complexes_subunits_WT_vs_KO.pdf", sep='')
                pdf(pdf.name, width=1.8, height=1.6)
                par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
                
                for(iindex in index)
                {
                    name = nuclear[iindex,20]
                    y0 = as.numeric(nuclear[iindex, c(1:16)]);
                    
                    averg = mean.err(y0, period=24)[1,]
                    err= mean.err(y0, period=24)[2,]

                    kk = grep('Cry.KO', colnames(nuclear))
                    y1 = as.numeric(nuclear[iindex, kk])
                    time = c(0:7)*3
                    lims = range(c(averg-err, averg+err, y1), na.rm=TRUE)
                    
                    plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=toupper(name), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
                    points(time, averg, col="darkblue", lwd=1.2, type='l')
                    points(6*c(0:3), y1, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
                    cex = 0.8
                    axis(1,at=6*c(0:4),cex.axis =cex)
                    axis(2,cex.axis =cex)
                    
                    box()
                    
                    abline(h=0,lty=2,lwd=1.5, col="darkgray")
        
                }
                dev.off()
            }
        }
        
        #### proteins invovled in cell cycle for WT and KO
        gg = read.delim('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/genes_cell_cycle.txt', sep='\t', header=FALSE)
        gg = unique(gg[,3])
        index = nuclear.names[match(gg, nuclear.names[,3]),1]
        index = index[which(nuclear$pval[index]<0.05)]
        
        pdf.name = paste("myplots/FIGURES/cell_markers_KO/Rhythmic_proteins_cell_cycle_WT_vs_KO.pdf", sep='')
        pdf(pdf.name, width=1.8, height=1.6)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        for(iindex in index)
        {
            name = nuclear[iindex,20]
            y0 = as.numeric(nuclear[iindex, c(1:16)]);
            
            averg = mean.err(y0, period=24)[1,]
            err= mean.err(y0, period=24)[2,]
            
            kk = grep('Cry.KO', colnames(nuclear))
            y1 = as.numeric(nuclear[iindex, kk])
            time = c(0:7)*3
            lims = range(c(averg-err, averg+err, y1), na.rm=TRUE)
            
            plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=toupper(name), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
            points(time, averg, col="darkblue", lwd=1.2, type='l')
            points(6*c(0:3), y1, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
            cex = 0.8
            axis(1,at=6*c(0:4),cex.axis =cex)
            axis(2,cex.axis =cex)
            
            box()
            
            abline(h=0,lty=2,lwd=1.5, col="darkgray")
            
        }
        dev.off()


        #### Targets of kinases CKD1, CKD4, CDK6 and Wee1
        load(file='Rdata/Phos_Nuclear_Mapping_all_v2.Rdata')
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/functions_nuclear.R')
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
        
        table.sx = read.table(file='Tables_DATA/Table_Nuclear_Prot_v3.txt', sep='\t', header=TRUE, as.is = c(2,3,10:19))
        nuclear = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/nuclear_proteins_L_H_log2_all_WT_KO_24h_12h_statistics.txt', sep='\t', header=TRUE, as.is=c(17:20))
        nuclear.names = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/Gene_names_Mapping_Nuclear_proteins.txt',header=TRUE, sep='\t')
        load(file='Rdata/Phos_all_with_seq_names_v2.Rdata')
        
        load(file='Rdata/kinase_motifs_curated_occurrence_matrix_v2.Rdata')
        
        load(file='Rdata/Kinase_candidates_elastic_net_alpha_0.02_curated.Rdata')
        load(file='Rdata/Kinase_candidates_phase_enrichment_p_0.05_curated.Rdata')
        load(file='Rdata/kinases_phosphatases_keep.Rdata')
        
        load(file='Tables_DATA/Kinases_Candidates_list_curated_targets.Rdata')
        
        mat.oc = mat.ocm
        kmotifs = kmotifs
        kinase.candidates = kinase.candidates
        infer = kinase.candidates[,-1]
        
        standardization.nona = function(x)
        {
            xx = (x-mean(x[which(!is.na(x)==TRUE)]))/sd(x[which(!is.na(x)==TRUE)])
            return(xx)
        }
        
        rsum = apply(mat.oc, 1, sum)
        
        mtfs = c("Cdk1_ST", "Cdk3_4_5_6_ST", "Csnk1g3_Rps6ka4_Wee1_ST")
        kk = match(mtfs, rownames(infer))
        
        ref = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Database_PhosphoNetworks/rawKSI.csv',sep='\t', header=TRUE)
        
        for(n.index in 1:length(kk))
        {
            n = kk[n.index];
            
            mm = match(rownames(infer)[n], colnames(mat.oc))
            nn = which((mat.oc[,mm])>0)
            targets = rownames(mat.oc)[nn]
            index = match(targets, phospho$seq.names)
            index = index[order(phospho$pval[index])]
            
            pdf.name = paste("myplots/FIGURES/cell_markers_KO/Targets_Kinases_", rownames(infer)[n],"_WT_vs_KO.pdf", sep = "")
            #pdf.name = paste("myplots/FIGURES/cell_markers_KO/", pc.name, "_Rhythmic_Complexes_subunits_WT_vs_KO.pdf", sep='')
            pdf(pdf.name, width=10, height=2.5)
            par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            for(iindex in index)
            {
                par(mfrow=c(1,3))
                #name = phospho$gene[iindex]
                name = phospho$seq.names[iindex];
                y0 = as.numeric(phospho[iindex, c(1:16)]);
                
                averg = mean.err(y0, period=24)[1,]
                err= mean.err(y0, period=24)[2,]
                
                ko.index = grep('Cry.KO', colnames(phospho))
                y1 = as.numeric(phospho[iindex, ko.index])
                time = c(0:7)*3
                lims = range(c(averg-err, averg+err, y1), na.rm=TRUE)
                
                plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=paste(toupper(name), ' \npval = ', signif(phospho$pval[iindex], d=3)), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
                points(time, averg, col="darkblue", lwd=1.2, type='l')
                points(6*c(0:3), y1, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
                cex = 0.8
                axis(1,at=6*c(0:4),cex.axis =cex)
                axis(2,cex.axis =cex)
                
                box()
                
                abline(h=0,lty=2,lwd=1.5, col="darkgray")
                
                
                gg = phospho.nuclear$mapping.gene[iindex];
                jj = which(nuclear$Gene.names==gg)
                jj= jj[which(nuclear$pval[jj]==min(nuclear$pval[jj]))]
                
                y00 = as.numeric(nuclear[jj, c(1:16)]);
                averg = mean.err(y00, period=24)[1,]
                err= mean.err(y00, period=24)[2,]
                
                ko.index = grep('Cry.KO', colnames(nuclear))
                y11 = as.numeric(nuclear[jj, ko.index])
                time = c(0:7)*3
                lims = range(c(averg-err, averg+err, y11), na.rm=TRUE)
                
                plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=paste(toupper(nuclear$Gene.names[jj])), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
                points(time, averg, col="darkblue", lwd=1.2, type='l')
                points(6*c(0:3), y11, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
                cex = 0.8
                axis(1,at=6*c(0:4),cex.axis =cex)
                axis(2,cex.axis =cex)
                
                box()
                
                abline(h=0,lty=2,lwd=1.5, col="darkgray")
                
                
                y0.ratio = y0/y00;
                averg = mean.err(y0.ratio, period=24)[1,]
                err= mean.err(y0.ratio, period=24)[2,]
                
                y1.ratio = y1/y11;
                time = c(0:7)*3
                lims = range(c(averg-err, averg+err, y1.ratio), na.rm=TRUE)
                
                plotCI(time, averg, err, scol="darkblue",lwd=1.2, cex=0.8, pch=16, col="darkblue", main=paste(toupper(nuclear$Gene.names[jj])), ylim=lims, ylab=NA, cex.lab=1.0, cex.main=0.8, xlab=NA, xlim=c(0,24), axes=FALSE)
                points(time, averg, col="darkblue", lwd=1.2, type='l')
                points(6*c(0:3), y1.ratio, type='b', pch=16, lwd=1.2, cex=0.8, col="darkred")
                cex = 0.8
                axis(1,at=6*c(0:4),cex.axis =cex)
                axis(2,cex.axis =cex)
                
                box()
                
                abline(h=0,lty=2,lwd=1.5, col="darkgray")
                
            }
            dev.off()

        }

        
        #### Targets of E2f1
        E2f1_targets = FALSE
        if(E2f1_targets)
        {
            load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/Elastic_TFs_Input_Pol2_CM_library_final.Rdata')
            motifs = as.matrix(mot_pol2)
            kk = grep('E2F', colnames(motifs))
            nb.targets = motifs[which(motifs[,kk]>0), kk]
            names.targets = rownames(motifs)[which(motifs[,kk]>0)]
            jj = match(names.targets, rownames(res.pol2.sel))
            jj = jj[which(res.pol2.sel[jj, 6]<0.01)]
            hist(res.pol2.sel[jj, 5])
            
            ensgene = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/mapping_ens_gene_trans.txt', sep='\t', header=FALSE)
            gene.mapping = read.table('/Users/jiwang/RNA_seq_Data/usful_mapping_files/genes_mapping.txt', sep=' ', header=FALSE)
            mm = match(names.targets, ensgene[,2])
            targets = ensgene[mm, 1]
            mm = match(targets, gene.mapping[,1])
            targets = gene.mapping[mm, 2]
            
            Filter.By.Chip = FALSE
            if(Filter.By.Chip)
            {
                
            }
            kk = match(targets, nuclear.names[,3])
            index = nuclear.names[kk, 1]
            index = index[which(!is.na(index)==TRUE)]
            mm = match(index, table.sx[,1])
            mm = mm[which(!is.na(mm)==TRUE)]
            mm = mm[which(table.sx$Localization[mm]=='Nucleus'|table.sx$Localization[mm]=='Nucleus/Cytoplasm')]
            
            pdf.name = paste("myplots/FIGURES/Targets_E2fs.pdf", sep = "")
            pdf(pdf.name, width=1.7, height=1.5)
            par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            #index = index[which(phospho$qv[index]<0.05)]
            
            mm.index = table.sx[mm, 1]
            test = as.matrix(nuclear[mm.index, c(1:16)])
            test = t(apply(test, 1, standardization.nona))
            
            time = seq(0, 45, by=3)
            
            lims = range(test, na.rm=TRUE)
            
            matplot(c(0:15)*3, t(test), type='n', col="darkblue",lwd=1.2,cex=0.6, main='E2Fs_motif', cex.main=1.2, ylim=lims, ylab=NA, xlab=NA, xlim=c(0,48), axes=FALSE)
            for(ii in 1:nrow(test))
            {
                points(time, test[ii, ], type='l', lwd=1.0, col='lightblue3')
            }
            
            ss = apply(test, 2, mean, na.rm=TRUE)
            points(time, ss, type='l', lwd=2.0, col='black')
            
            axis(1,at=seq(0, 48, by=12),cex.axis =1.0)
            #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
            lims = signif(lims, d=1)
            if(gene=='Hdac3') {
                axis(2,at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2),las=1,cex.axis = 1.0)
            }else{
                by = signif((lims[2]-lims[1])/5,d=0)
                axis(2,at = seq(lims[1], lims[2], by=by),las=1,cex.axis = 1.0)
            }
            box()
            abline(h=0,lty=2,lwd=1.2, col="darkgray")
            #if(gene=='Rorc')
            #legend(0,(lims[2]+0.5),c('phospho','nuclear'),pch=c(16,16),lty=c(1,1),col=c("darkgreen", "darkblue"), bty='n')
            
            dev.off()
        }
        
    }
    ########
    #### Examples of rhythmic proteins involved in cell cycle
    ########
    hmr = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/DATA_WB_2.csv', sep=',', header=TRUE)
    colnames(hmr) = toupper(colnames(hmr))
    source('functions_nuclear.R')
    
    require('plotrix')
    examples2 = c('Mcm2', 'Mcm3', 'P.H3.Ser10', 'Rb', 'P.rb.ser807.811')
    
    match(examples2, nuclear.names[,3])
    match(examples2, phospho.names[,2])
    
    for(gene in examples2)
    {
        pdf.name = paste("myplots/FIGURES/Fig_examples_4cell_cycle_",gene,".pdf", sep = "")
        pdf(pdf.name, width=1.7, height=1.7)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        kk = grep(toupper(gene), colnames(hmr))
        if(gene=='P.H3.Ser10')
        kk = c(kk, NA, NA)
        y1 = c(as.numeric(unlist(hmr[c(1:16), kk[1]])), as.numeric(unlist(hmr[c(1:16), kk[2]])))
        y1 = y1/mean(y1, na.rm=TRUE)
        #y1 = rbind(y1[c(1:8)], y1[c(9:16)], y1[17:24], y1[25:32])
        y2 = c(as.numeric(unlist(hmr[c(1:16), kk[3]])), as.numeric(unlist(hmr[c(1:16), kk[4]])))
        y2 = y2/mean(y2, na.rm=TRUE)
        
        averg1 = mean.err(y1, period=24)[1,]
        err1 = mean.err(y1, period=24)[2,]
        averg2 = mean.err(y2, period=24)[1,]
        err2 = mean.err(y2, period=24)[2,]
        
        time = c(0:7)*3
        lims = range(c(averg2-err2, averg2+err2, averg1-err1, averg1+err1), na.rm=TRUE)
        
        if(gene=='Mcm2'|gene=='Mcm3')
        {
            jj = match(gene, nuclear.names[,3])
            y0 = as.numeric(nuclear[nuclear.names[jj, 1], c(1:16)])
            y0 = 2^y0
            y0 = y0/mean(y0, na.rm=TRUE)
            averg0 = mean.err(y0, period=24)[1,]
            err0 = mean.err(y0, period=24)[2,]
            lims = range(c(averg0-err0, averg0+err0, averg1-err1, averg1+err1), na.rm=TRUE)
        }
        
        #lims = range(c(y0[which(!is.na(y0)==TRUE)]))
        plot(time, averg1, type='l', ylim=lims, col='red', lwd=1.2, pch=16, main=NA, axes=FALSE, xlab=NA, ylab=NA)
        points(time, averg1, type='p', cex=0.7, col='red', pch=16)
        arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3, col='red', lwd=1.2, cex=0.7)
        
        points(time, averg2, type='l', lwd=1.2, col='black', pch=16)
        points(time, averg2, type='p', cex=0.7, col='black', pch=16)
        arrows(time, averg2-err2, time, averg2+err2, length=0.05, angle=90, code=3, col='black', lwd=1.2, cex=0.7)
        
        if(gene=='Rb'|gene=='P.rb.ser807.811') axis(1,at=6*c(0:4),cex.axis =1)
        #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
        lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
        by = 0.1
        #print(gene)
        #print(lims)
        if(gene=='Rb') axis(2,at = c(0.8, 1.0, 1.2, 1.4),las=1,cex.axis = 1)
        if(gene=='Mcm2') axis(2,at = c(0.5, 1.0, 1.5),las=1,cex.axis = 1)
        if(gene=='Mcm3') axis(2,at = c(0.5, 1.0, 1.5),las=1,cex.axis = 1)
        if(gene=='P.H3.Ser10') axis(2,at = c(0.6, 1.0, 1.4),las=1,cex.axis = 1)
        if(gene=='P.rb.ser807.811') axis(2,at = c(0.5, 1.0, 1.5, 2.0),las=1,cex.axis = 1)
        #axis(2,at = seq(lims[1], lims[2],by=0.3),las=1,cex.axis = 1)
        
        #if(gene=='Arntl') axis(2,at = c(-0.6, -0.3, 0.0, 0.3),las=1,cex.axis = 1.2)
        
        box()
        abline(h=1,lty=2,lwd=1.0, col="darkgray")
        
        if(gene=='Mcm2'|gene=='Mcm3')
        {
            points(time, averg0, type='l', lwd=1.2, col='darkblue', pch=15)
            points(time, averg0, type='p', pch=15, col='darkblue', cex=0.7)
            arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col='darkblue', lwd=1.2)
            #legend(13, (lims[2]), c('MS.NE','WB.NE', 'WB.Cyto'),pch=c(16,16, 15),lty=c(1,1, 1),col=c("darkblue","red", "black"), cex=0.7, pt.lwd=1.5, pt.cex=1.0, bty='n')
            
        }
        
        dev.off()
    }
    
    
    ########
    #### % of liver cells in cell cycle marked by Ki-67
    ########
    #load(file='Rdata/Images_results_KI67_marker.Rdata')
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    load(file='Rdata/Images_results_KI67_marker_all_v4.Rdata')
    source('functions_nuclear.R')
    positive = keep
    
    index = c(1:16, 19:34)
    kk = match(index, keep[,1])
    #positive[,4] = positive[,4]
    averg0 = mean.err(positive[kk, 5])[1,]
    err0 = mean.err(positive[kk, 5])[2,]
    
    lims = range(c(averg0+err0, averg0-err0));
    
    ttime = c(0:7)*3
    tt = seq(0, 24, length.out=100)
    
    pdf('myplots/FIGURES/Ki-67_positive_nuclei.pdf', width=2.2, height=1.0)
    #pdf('myplots/FIGURES/Cell_DNA_contents_8N.pdf', width=2.2, height=1.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2.0,3,0.5,0.8)+0.1, tcl = -0.3)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,2,0.8,0.8)+0.1, tcl = -0.3)
    #par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2.0,2.0,0.5,0.8)+0.1, tcl = -0.3)
    
    pch=16
    col = 'plum1'
    #col2 = 'gray'
    lims = c(0, 10)
    lwd = 1.2
    cex = 1.0
    cex2 = 0.6
    plot(ttime, averg0, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd, cex=cex)
    res = f24_R2_alt2(as.numeric(positive[index,4]), t=c(0:31)*3)
    #yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    #apoints(tt, yy, type='l', col=col, lwd=lwd)
    
    #points(ttime, mean.err(kkeep[,11])[1,], type='p', pch=pch, col=col2, cex=cex2, lwd=lwd)
    #averg = mean.err(kkeep[,11])[1,]
    #err = mean.err(kkeep[,11])[2,]
    #arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col2, lty=1, pch=pch, lwd=lwd, cex=cex)
    abline(h=mean(positive[kk,5]), col=col)
    axis(2, at= seq(0,10, by=5), las=1,cex.axis = 1.0)
    axis(1,at=6*c(0:4),cex.axis =1.0)
    box()
    
    dev.off()
    
    Size.KO = FALSE
    if(Size.KO)
    {
        #### median of nuclei sizes
        averg0 = mean.err(positive[kk, 10])[1,]
        err0 = mean.err(positive[kk, 10])[2,]
        ttime = c(1:8)*3
        
        averg0 = averg0[c(2:8, 1)]
        err0 = err0[c(2:8, 1)]
        #ttime = ttime[c(2:8, 1)]
        
        #tt = seq(0, 24, length.out=100)
        
        pdf('myplots/FIGURES/Ki-67_positive_nuclei_median_size.pdf', width=2., height=1.7)
        par(cex = 0.7, las = 1, mgp = c(0.5,0.5,0), mar = c(2,2,0.8,0.8)+0.1, tcl = -0.3)
        
        pch=21
        col = 'darkmagenta'
        #col2 = 'gray'
        #lims = c(0, 10)
        lims = range(c(averg0+err0, averg0-err0));
        lwd = 1.2
        cex = 1.0
        cex2 = 0.6
        plot(ttime, averg0, type='l', ylim=lims, xlim=c(3,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
        points(ttime, averg0, type='p', col=col, pch=pch)
        arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd, cex=cex)
        res = f24_R2_alt2(as.numeric(positive[index,4]), t=c(0:31)*3)
        #yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
        #apoints(tt, yy, type='l', col=col, lwd=lwd)
        
        #points(ttime, mean.err(kkeep[,11])[1,], type='p', pch=pch, col=col2, cex=cex2, lwd=lwd)
        #averg = mean.err(kkeep[,11])[1,]
        #err = mean.err(kkeep[,11])[2,]
        #arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col2, lty=1, pch=pch, lwd=lwd, cex=cex)
        #abline(h=mean(positive[kk,5]))
        axis(2, at= seq(180,240, by=20), las=1,cex.axis = 1.0)
        axis(1,at=ttime,cex.axis =1.0)
        box()
        
        dev.off()
        
        #### WT vs Cry DKO
        source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
        load(file='Rdata/Images_results_KI67_marker_all_v4.Rdata')
        positive = keep
        positive = positive[order(positive[,1]), ]
        #index = c(1:16, 19:34)
        x1 = positive[c(4,12,22,30), 5]
        x2 = positive[c(8,16,26,34), 5]
        x11 = positive[c(17, 35), 5]
        x22 = positive[c(18, 36), 5]
        
        t.test(x1, x11)
        t.test(x2, x22)
        t.test(x1, x2)
        
        ym = c(mean(x1), mean(x11), NA, mean(x2), mean(x22))
        dy = c(sd(x1)/sqrt(4), sd(x11)/sqrt(2), NA, sd(x2)/sqrt(4), sd(x22)/sqrt(2))
        cols = c('gray70', 'gray30', 'white', 'gray70', 'gray30')
        
        pdf('myplots/FIGURES/Ki-67_positive_percentages_WT_KO.pdf', width=2.0, height=1.4)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,3,0.5,0.8)+0.1, tcl = -0.3)
        
        pch = 16
        barx = barplot(ym, col=cols, ylim=c(0, 10), axes=FALSE, w=0.2, lwd=0.01)
        arrows(barx, ym-dy, barx, ym+dy, length=0.04, angle=90, code=3, col='black', lty=1, pch=pch,lwd=1.2)
        axis(2, at= seq(0,6,by=2), las=1,cex.axis = 1.0)
        
        dev.off()
    }
    
    #############
    ## DNA contents by FACS and estimations from nuclei sizes
    #############
    add.KO = FALSE
    dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
    #load(file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32.Rdata')
    source('functions_nuclear.R')
    
    ref2 = dna[1, c(2:33)]
    ref4 = dna[2, c(2:33)]
    ref8 = dna[4, c(2:33)]
    averg0 = mean.err(ref2)[1,]
    err0 = mean.err(ref2)[2,]
    averg1 = mean.err(ref4)[1,]
    err1 = mean.err(ref4)[2,]
    averg2 = mean.err(ref8)[1,]
    err2 = mean.err(ref8)[2,]
    lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1, averg2+err2, averg2-err2));
    
    ttime = c(0:7)*3
    tt = seq(0, 24, length.out=100)
    
    pdf('myplots/FIGURES/Cell_DNA_contents_2N.pdf', width=2.2, height=0.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    pch=16
    col = 'black'
    col2 = 'gray'
    lims = c(42, 70)
    lwd = 1.2
    cex = 1.0
    cex2 = 0.6
    plot(ttime, averg0, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd, cex=cex)
    res = f24_R2_alt2(as.numeric(ref2), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=lwd)
    if(add.KO)
    {
        ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[1, c(2:3)]), as.numeric(ko[1, c(4:5)]))
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
    }
    
    #points(ttime, mean.err(kkeep[,11])[1,], type='p', pch=pch, col=col2, cex=cex2, lwd=lwd)
    #averg = mean.err(kkeep[,11])[1,]
    #err = mean.err(kkeep[,11])[2,]
    #arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col2, lty=1, pch=pch, lwd=lwd, cex=cex)
    
    axis(2, at= seq(50, 70, 10), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    pdf('myplots/FIGURES/Cell_DNA_contents_4N.pdf', width=2.2, height=0.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = c(30, 50)
    
    plot(ttime, averg1, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=1.0,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg1-err1, ttime, averg1+err1, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch, cex=cex)
    res = f24_R2_alt2(as.numeric(ref4), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=lwd)
    if(add.KO)
    {
        ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[2, c(2:3)]), as.numeric(ko[2, c(4:5)]))
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
    }
    
    #points(ttime, mean.err(kkeep[,12])[1,], type='p', pch=pch, col=col2, cex=cex2)
    #averg = mean.err(kkeep[,12])[1,]
    #err = mean.err(kkeep[,12])[2,]
    #arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col2, lty=1, pch=pch, lwd=lwd, cex=cex2)
    
    axis(2, at= seq(30,50, 10), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    
    pdf('myplots/FIGURES/Cell_DNA_contents_8N.pdf', width=2.2, height=1.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2.0,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = c(0, 10)
    plot(ttime, averg2, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=1.0,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg2-err2, ttime, averg2+err2, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    res = f24_R2_alt2(as.numeric(ref8), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=lwd)
    if(add.KO)
    {
        ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[4, c(2:3)]), as.numeric(ko[4, c(4:5)]))
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
    }

    #points(ttime, mean.err(kkeep[,13])[1,], type='p', pch=pch, col=col2, cex=cex2)
    #averg = mean.err(kkeep[,13])[1,]
    #err = mean.err(kkeep[,13])[2,]
    #arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col2, lty=1, pch=pch, lwd=lwd, cex=cex2)
    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    axis(2, at= c(0,5, 10), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    #############
    ####### Results of image analysis
    #############
    load(file='Rdata/Results_images_WT_vf_calibrated_FACS_outliers_1_22_29_32.Rdata')
    add.KO = FALSE
    #load(file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32.Rdata')
    
    ##### test nb of cells and size of images
    plot(c(0:31)*3, keep[,1]*keep[,2]/max(keep[,1]*keep[,2]), type='b', ylim=c(0, 1));
    points(c(0:31)*3, keep[,3]/max(keep[,3]), type='b', col='blue')
    abline(v=c(1:3)*24, col='red')
    abline(v=(c(0:3)*24+9), col='gray')
    
    source('functions_nuclear.R')
    
    remv = c(1, 22, 29, 32)
    kkeep = keep
    kkeep[remv, ] = NA
    
    Table.Sup = FALSE
    if(Table.Sup)
    {
        ### motifs
        xx = dna
        xx[,1] = c('2N', '4N', '6N', '8N')
        colnames(xx)[1] = 'DNA.content'
        colnames(xx)[2:9] = paste('ZT', c(0:7)*3, '.rep1', sep='')
        colnames(xx)[10:17] = paste('ZT', c(0:7)*3, '.rep2', sep='')
        colnames(xx)[18:25] = paste('ZT', c(0:7)*3, '.rep3', sep='')
        colnames(xx)[26:33] = paste('ZT', c(0:7)*3, '.rep4', sep='')
        colnames(xx)[34] = 'nb.timepoints'
        xx = xx[, -40]
        #colnames(xx) = c('motifs', 'phase.motifs', 'amp.motifs', 'TFs.associated', 'phase.TFs', 'pval.TFs', 'qv.TFs')
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S7_DNA_contents_FACS.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
     
        xx = kkeep
        xx = t(xx)
        xx = data.frame(xx, stringsAsFactors=FALSE)
        colnames(xx)[c(1:8)] = paste('ZT', c(0:7)*3, '.rep.1', sep='')
        colnames(xx)[c(9:16)] = paste('ZT', c(0:7)*3, '.rep.2', sep='')
        colnames(xx)[c(17:24)] = paste('ZT', c(0:7)*3, '.rep.3', sep='')
        colnames(xx)[c(25:32)] = paste('ZT', c(0:7)*3, '.rep.4', sep='')
        
        jj = c(3, 4, 5, 6, 14, 15, 17, 16, 18)
        xx = xx[jj,]
        rownames(xx)[c(3:9)] = c('cell.area.median.mono',  'cell.area.median.bi', 'cell.percent.2N.mono', 'cell.percent.4N.mono', 'cell.percent.4N.bi', 'cell.percent.8N.mono', 'cell.percent.8N.bi')
        xx[c(3:4), ] =  xx[c(3:4), ]*scale
        
        write.table(xx, file='/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Paper/Supplemental_Tables/Table_S7_summary_image_analysis.txt', sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)
    }
    
    ref2 = dna[1, c(2:33)]/100
    ref4 = dna[2, c(2:33)]/100
    ref8 = dna[4, c(2:33)]/100
    
    ########
    #### Cell size
    ########
    pdf('myplots/FIGURES/Cell_Sizes_nbs.pdf', width=2.2, height=1.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    source('functions_nuclear.R')
    #averg0 = mean.err(keep$cell.size.all.mean, period=48)[1,]
    #err0 = mean.err(keep$cell.size.all.mean, period=48)[2,]
    scale = 0.34*0.34
    averg0 = mean.err(keep$cell.size.mono.mean*scale, period=24)[1,]
    err0 = mean.err(keep$cell.size.mono.mean*scale, period=24)[2,]
    averg1 = mean.err(keep$cell.size.binu.mean*scale, period=24)[1,]
    err1 = mean.err(keep$cell.size.binu.mean*scale, period=24)[2,]
    
    lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1));
    lims = c(180, 500)
    pch=16
    
    time = c(0:7)*3
    xlim = c(0, 24)
    col = 'blue'
    plot(time, averg0, type='l', ylim=lims, xlim=xlim, col=col, lwd=1.2,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    points(time, averg0, type='p', col=col, cex=1., pch=pch)
    arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.2)
    col='red'
    pch = 16 ;
    points(time, averg1, type='l', col=col, lwd=1.2, lty=1)
    points(time, averg1, type='p', col=col, lwd=1.2, lty=1, pch=pch, cex=1.0)
    arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3,  lty=1, col=col, pch=pch, lwd=1.2)
    #points(time, averg2, type='l', col=col, lwd=1.5,  lty=4, pch=pch)
    #arrows(time, averg2-err2, time, averg2+err2, length=0.05, angle=90, code=3,  lty=1,col=col, pch=pch,lwd=1.5)
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 5]), as.numeric(ko[c(3:4), 5]))*scale
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 6]), as.numeric(ko[c(3:4), 6]))*scale
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)

    }

    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2, at= seq(200, 500, by=100), las=1,cex.axis = 1.0)
    box()
    #legend(7, 4500, c('mononu', 'binu'), pch=c(1, 0), lty=c(1,2),col=c('blue', 'red'), lwd=1.5, cex=0.9, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()
    
    ######
    ### Glycogen in WT and KO
    ######
    glyco = read.table('/Users/jiwang/Dropbox/GachonProt/Nuclear_Prot/Glycogen_WT_CRYKO/GLYCOGEN.txt', sep='\t', header=FALSE)
    source('functions_nuclear.R')
    
    averg0 = mean.err(glyco[3, c(2:19)], period=24, interval=4)[1,]
    err0 = mean.err(glyco[3, c(2:19)], period=24, interval=4)[2,]
    averg1 = mean.err(glyco[2, c(2:19)], period=24, interval=4)[1,]
    err1 = mean.err(glyco[2, c(2:19)], period=24, interval=4)[2,]
    
    lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1));
    lims = c(0, 100)
    pch=16
    time = c(0:5)*4
    xlim = c(0, 24)
    
    pdf('myplots/FIGURES/glycogen_variation_WT_KO.pdf', width=2.2, height=1.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    col = 'black'
    pch = 16
    plot(time, averg0, type='l', ylim=lims, xlim=xlim, col=col, lwd=1.2,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    points(time, averg0, type='p', col=col, cex=1., pch=pch)
    arrows(time, averg0-err0, time, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.2)
    col='gray70'
    pch = 21;
    #points(time, averg1, type='l', col=col, lwd=1.2, lty=1)
    #points(time, averg1, type='p', col=col, lwd=1.2, lty=1, pch=pch, cex=1.0)
    #arrows(time, averg1-err1, time, averg1+err1, length=0.05, angle=90, code=3,  lty=1, col=col, pch=pch, lwd=1.2)
    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2, at= seq(0, 100, by=50), las=1,cex.axis = 1.0)
    box()
    #legend(7, 4500, c('mononu', 'binu'), pch=c(1, 0), lty=c(1,2),col=c('blue', 'red'), lwd=1.5, cex=0.9, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()

    
    ##################
    ###### Binucleated cell percentages
    ##################
    #filename = paste('Rdata/Results_images_KO_v2.Rdata', sep='')
    #load(file=filename)
    pdf('myplots/FIGURES/Cell_percent_mono.pdf', width=2., height=0.8)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3.0,0.5,0.8)+0.1, tcl = -0.3)
    
    averg = mean.err((1-keep$percent.binucleated), period=24)[1,]
    err = mean.err((1-keep$percent.binucleated), period=24)[2,]
    
    col = "blue"
    #col = 'black'
    #lims = range(c(averg+err, averg-err))
    pch = 16
    lims = c(0.8, 0.9)
    plot(time, averg, type='l', ylim=lims, xlim=c(0, 21), col=col, lwd=1.2,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    points(time, averg, type='p', col=col, cex=1.0, pch=pch)
    arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.2)
    axis(2, at= c(0.8, 0.9), las=1,cex.axis = 1.0)
    #legend(15, 0.91, c('mono'), pch=c(1), lty=c(1),col=c('blue'), lwd=1.2, cex=1.0, pt.lwd=1.2, pt.cex=1.0, bty='n')
    box()
    
    dev.off()
    
    
    pdf('myplots/FIGURES/Cell_percent_bino.pdf', width=2.2, height=1.4)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    averg = mean.err(keep$percent.binucleated*100, period=24)[1,]
    err = mean.err(keep$percent.binucleated*100, period=24)[2,]
    
    #averg1 = c(mean(keep.ko$percent.binucleated[c(1,3)]), mean(keep.ko$percent.binucleated[c(2,4)]))
    #err1 = c(sd(keep.ko$percent.binucleated[c(1,3)]), sd(keep.ko$percent.binucleated[c(2,4)]))/sqrt(2)
    
    col = "red"
    #col = 'black'
    #lims = range(c(averg+err, averg-err))
    pch = 16
    lims = c(10, 20)
    plot(time, averg, type='n', ylim=lims, xlim=c(0, 24), col=col, lwd=1.2,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    points(time, averg, type='p', col=col, cex=1.0, pch=16)
    arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.2)
    
    res = f24_R2_alt2(as.numeric(keep$percent.binucleated*100), t=c(0:31)*3)
    tt = seq(0, 21, by=0.1)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=lwd)
    
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 4]), as.numeric(ko[c(3:4), 4]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }
    #points(c(9,21), averg1, type='p', col='darkgray', cex=1.0)
    #arrows(c(9,21), averg1-err1, c(9,21), averg1+err1, length=0.05, angle=90, code=3, col='darkgray', lty=1, pch=pch,lwd=1.2)
    axis(2, at= c(10, 15, 20), las=1,cex.axis = 1.0)
    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    box()
    #legend(15, 0.21, c('binu'), pch=pch, lty=c(2),col=c('red'), lwd=1.2, cex=0.9, pt.lwd=1.2, pt.cex=1.0, bty='n')
    
    dev.off()

    #############
    ####### cell density
    #############
    pdf('myplots/FIGURES/Cell_density.pdf', width=2.2, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    ss = keep$image.width*keep$image.length
    ss = ss/min(ss)
    averg = mean.err(keep$nb.cell.detected/ss ,period=24)[1,]
    err = mean.err(keep$nb.cell.detected/ss, period=24)[2,]
    
    #col = "cadetblue4"
    col = 'black'
    lims = range(c(averg+err, averg-err))
    plot(time, averg, type='l', ylim=lims, xlim=c(0,24), col=col, lwd=1.5,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA)
    arrows(time, averg-err, time, averg+err, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.5)
    axis(1,at=3*c(0:8),cex.axis =1.0)
    axis(2, at= seq(6000, 12000, by=3000), las=0,cex.axis = 1.0)
    box()
    #legend('topright', legend = c('Mauvoisin','Robles'), lty=c(1,1), cex=0.7,col = c('blue', 'green'), border = NA, bty = 'n')
    dev.off()
    
    #############
    ### nucleus size for different DNA centents
    #############
    ttime = c(0:7)*3
    #tt = seq(0, 24, length.out=100)
    averg0 = mean.err(kkeep[,8]*scale)[1,]
    err0 = mean.err(kkeep[,8]*scale)[2,]
    
    averg1 = mean.err(kkeep[,9]*scale)[1,]
    err1 = mean.err(kkeep[,9]*scale)[2,]
    
    averg2 = mean.err(kkeep[,10]*scale)[1,]
    err2 = mean.err(kkeep[,10]*scale)[2,]
    
    lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1, averg2+err2, averg2-err2));
    
    pdf('myplots/FIGURES/Nucleus_sizes.pdf', width=2.2, height=1.6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    #col = 'black'
    col = 'darkmagenta'
    cex = 1.0
    lims = c(20, 110)
    pch=16
    lwd = 1.2
    #lims = c(0.5, 0.7)
    plot(ttime, averg0, type='l', ylim=lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    points(ttime, averg0, type='p', col=col, cex=cex, pch=pch)
    arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    #res = f24_R2_alt2(as.numeric(ref2), t=c(0:31)*3)
    #yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    #points(tt, yy, type='l', col=col, lwd=1.)
    
    pch = 16
    points(ttime, averg1, type='l', pch=pch, col=col, cex=cex)
    points(ttime, averg1, type='p', col=col, cex=cex, pch=pch)
    arrows(ttime, averg1-err1, ttime, averg1+err1, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch, lwd=lwd)
    
    pch = 16
    points(ttime, averg2, type='l', pch=pch, col=col, cex=cex)
    points(ttime, averg2, type='p', col=col, cex=cex, pch=pch)
    arrows(ttime, averg2-err2, ttime, averg2+err2, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch, lwd=lwd)
    
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 8]), as.numeric(ko[c(3:4), 8]))*scale
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 9]), as.numeric(ko[c(3:4), 9]))*scale
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
        ko = keep.ko
        ko = rbind(as.numeric(ko[c(1:2), 10]), as.numeric(ko[c(3:4), 10]))*scale
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
    }

    
    axis(1,at=6*c(0:4),cex.axis =1.0)
    axis(2, at= seq(20,100, by=20), las=1,cex.axis = 1.0)
    
    box()
    
    #legend('topleft', c('n.2N', 'n.4N.n', 'n.8N'), pch=c(1,0,2), lty=c(1,1,1),col=col, lwd=1.5, cex=0.9, pt.lwd=1.5, pt.cex=1.0, bty='n')
    
    dev.off()
    
    #############
    ### percentages of cells with polyploidy
    ############
    
    pdf('myplots/FIGURES/Cell_populations_1.pdf', width=2.2, height=0.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    #tt = seq(0, 24, length.out=100)
    averg = mean.err(kkeep[,14]*100)[1,]
    err = mean.err(kkeep[,14]*100)[2,]
    ttime = c(0:7)*3
    tt = seq(0, 24, by=0.1)
    
    pch=16
    col = 'lightslateblue'
    lwd = 1.2
    cex= 1.0
    
    lims = c(0.35, 0.6)*100
    
    plot(ttime, averg, type='p', ylim = lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    
    res = f24_R2_alt2(as.numeric(kkeep[,14]*100), t=c(0:31)*3)
    yy = res[2]*(1+0.0*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.2)
    print(res)
    
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 14]), as.numeric(ko[c(3:4), 14]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }
    
    axis(2, at= seq(0, 60, by=10), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    
    pdf('myplots/FIGURES/Cell_populations_2.pdf', width=2.2, height=0.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = c(28, 40)
    averg = mean.err(kkeep[,15]*100)[1,]
    err = mean.err(kkeep[,15]*100)[2,]
    
    plot(ttime, averg, type='p', ylim = lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    
    res = f24_R2_alt2(as.numeric(kkeep[,15]*100), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.5)
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 15]), as.numeric(ko[c(3:4), 15]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }
    print(res)
    axis(2, at= seq(0, 60, by=10), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    pdf('myplots/FIGURES/Cell_populations_3.pdf', width=2.2, height=0.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = c(5, 15)

    averg = mean.err(kkeep[,17]*100)[1,]
    err = mean.err(kkeep[,17]*100)[2,]
    plot(ttime, averg, type='p', ylim = lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    
    res = f24_R2_alt2(as.numeric(kkeep[,17]*100), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.5)
    print(res)
    
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 17]), as.numeric(ko[c(3:4), 17]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }

    #axis(1,at=3*c(0:6),cex.axis =1.0)
    axis(2, at= seq(5, 15, by=5), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    pdf('myplots/FIGURES/Cell_populations_4.pdf', width=2.2, height=0.7)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(0.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    #tt = seq(0, 24, length.out=100)
    averg = mean.err(kkeep[,16]*100)[1,]
    err = mean.err(kkeep[,16]*100)[2,]
    
    lims = c(0, 10)
    pch=16
    col = 'lightslateblue'
    plot(ttime, averg, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    
    res = f24_R2_alt2(as.numeric(kkeep[,16]*100), t=c(0:31)*3)
    yy = res[2]*(1+0.0*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.5)
    print(res)
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 16]), as.numeric(ko[c(3:4), 16]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }

    #axis(1,at=3*c(0:7),cex.axis =1.0)
    axis(2, at= c(0, 10), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    pdf('myplots/FIGURES/Cell_populations_5.pdf', width=2.2, height=0.85)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.5,3,0.5,0.8)+0.1, tcl = -0.3)
    
    #tt = seq(0, 24, length.out=100)
    
    averg = mean.err(kkeep[,18]*100)[1,]
    err = mean.err(kkeep[,18]*100)[2,]
    
    lims = c(0, 10)
    pch=16
    col = 'lightslateblue'
    plot(ttime, averg, type='p', ylim=lims, xlim=c(0,24), col=col, lwd=lwd,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=cex)
    arrows(ttime, averg-err, ttime, averg+err, length=0.04, angle=90, code=3, col=col, lty=1, pch=pch,lwd=lwd)
    
    res = f24_R2_alt2(as.numeric(kkeep[,18]*100), t=c(0:31)*3)
    yy = res[2]*(1+0.0*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.5)
    print(res)
    if(add.KO)
    {
        load(file='Rdata/Results_images_KO_calibrated_FACS.Rdata')
        ko = keep.ko
        #ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
        ko = rbind(as.numeric(ko[c(1:2), 18]), as.numeric(ko[c(3:4), 18]))*100
        
        averg.ko = apply(ko, 1, mean)
        err.ko = apply(ko, 1, sd)/sqrt(2)
        points(c(9, 21), averg.ko, type='p', pch=21, col='gray', cex=cex)
        arrows(c(9, 21), averg.ko-err.ko, c(9,21), averg.ko+err.ko, length=0.04, angle=90, code=3, col=col2, lty=1, pch=16, lwd=lwd, cex=cex)
        
    }

    axis(1,at=6*c(0:4),cex.axis =1.0)
    axis(2, at= c(0, 10), las=1,cex.axis = 1.0)
    box()
    dev.off()

    
    ##########
    ##### Compare WT and CRYDKO
    ##########
    
    #### compare DNA contents FACS
    load(file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32.Rdata')
    facs.ko = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
    index = c(4,12,20,28)+1
    ref21 = as.numeric(dna[1, index])
    ref41 = as.numeric(dna[2, index])
    ref81 = as.numeric(dna[4, index])
    ko21 = as.numeric(facs.ko[1,c(2:3)])
    ko41 = as.numeric(facs.ko[2,c(2:3)])
    ko81 = as.numeric(facs.ko[4,c(2:3)])
    
    index = c(8,16,24,32)+1
    ref22 = as.numeric(dna[1, index])
    ref42 = as.numeric(dna[2, index])
    ref82 = as.numeric(dna[4, index])
    ko22 = as.numeric(facs.ko[1,c(4:5)])
    ko42 = as.numeric(facs.ko[2,c(4:5)])
    ko82 = as.numeric(facs.ko[4,c(4:5)])
    
    t.test(ref21, ko21)
    t.test(ref22, ko22)
    #t.test(ref21, )
    t.test(ko21, ko22)
    wilcox.test(ko21, ko22)
    t.test(ref21, ref22)
    t.test(ref21[1:2], ref22[1:2])
    wilcox.test(ref21, ref22)
    wilcox.test(ref21, ref22, paires=TRUE, alternative = "two.sided")
    
    t.test(ref41, ko41)
    t.test(ref42, ko42)
    t.test(ko41, ko42)
    
    t.test(ref21, ref22)
    ym = c(
    
    mean(ref21), mean(as.numeric(facs.ko[1,c(2:3)])/100),
    mean(ref41), mean(as.numeric(facs.ko[2,c(2:3)])/100),
    mean(ref81), mean(as.numeric(facs.ko[4,c(2:3)])/100),
    NA,
    mean(ref22), mean(as.numeric(facs.ko[1,c(4:5)])/100),
    mean(ref42), mean(as.numeric(facs.ko[2,c(4:5)])/100),
    mean(ref82), mean(as.numeric(facs.ko[4,c(4:5)])/100)
    
    )
    
    dy = c(
    
    sd(ref21)/sqrt(4), sd(as.numeric(facs.ko[1,c(2:3)])/100)/sqrt(2),
    sd(ref41)/sqrt(4), sd(as.numeric(facs.ko[2,c(2:3)])/100)/sqrt(2),
    sd(ref81)/sqrt(4), sd(as.numeric(facs.ko[4,c(2:3)])/100)/sqrt(2),
    NA,
    sd(ref22)/sqrt(4), sd(as.numeric(facs.ko[1,c(4:5)])/100)/sqrt(2),
    sd(ref42)/sqrt(4), sd(as.numeric(facs.ko[2,c(4:5)])/100)/sqrt(2),
    sd(ref82)/sqrt(4), sd(as.numeric(facs.ko[4,c(4:5)])/100)/sqrt(2)

    )
    cols = c('gray70', 'gray30', 'gray70', 'gray30', 'gray70', 'gray30', 'white', 'gray70', 'gray30','gray70', 'gray30', 'gray70', 'gray30')
    
    pdf('myplots/FIGURES/DNA_Contents_FACS_WT_KO.pdf', width=2.5, height=1.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = range(c((ym-dy), (ym+dy)), na.rm=TRUE)
    lims[1] = 0
    barx = barplot(ym, col=cols, ylim=lims, axes=FALSE, w=0.2, lwd=0.01)
    arrows(barx, ym-dy, barx, ym+dy, length=0.04, angle=90, code=3, col='black', lty=1, pch=pch,lwd=1.2)
    
    axis(2, at= seq(0,0.6,by=0.2), las=1,cex.axis = 1.0)
    
    dev.off()
    
    #### compare binu cell percentagess
    filename = paste('Rdata/Results_images_KO_v2.Rdata', sep='')
    load(file=filename)
    load(file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32.Rdata')
    
    source('functions_images.R')
    source('functions_nuclear.R')
    
    x1 = keep$percent.binucleated[c(4,12,20,28)]
    x11 = keep.ko$percent.binucleated[c(1, 3)]
    x2 = keep$percent.binucleated[c(8,16,24,32)]
    x22 = keep.ko$percent.binucleated[c(2, 4)]
    t.test(x1, x11)
    t.test(x2, x22)
    t.test(x1, x2)
    
    
    ym = c(mean(keep$percent.binucleated[c(4,12,20,28)]), mean(keep.ko$percent.binucleated[c(1, 3)]), NA, mean(keep$percent.binucleated[c(8,16,24,32)]), mean(keep.ko$percent.binucleated[c(2, 4)]))
    dy = c(sd(keep$percent.binucleated[c(4,12,20,28)])/sqrt(4), sd(keep.ko$percent.binucleated[c(1, 3)])/sqrt(2), NA, sd(keep$percent.binucleated[c(8,16,24,32)])/sqrt(4), sd(keep.ko$percent.binucleated[c(2, 4)])/sqrt(2))
    cols = c('gray70', 'gray30', 'white', 'gray70', 'gray30')
    
    pdf('myplots/FIGURES/Binu_percents_WT_KO.pdf', width=1.6, height=1.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,3,0.5,0.8)+0.1, tcl = -0.3)
    
    pch = 16
    barx = barplot(ym, col=cols, ylim=c(0, 0.2), axes=FALSE, w=0.2, lwd=0.01)
    arrows(barx, ym-dy, barx, ym+dy, length=0.04, angle=90, code=3, col='black', lty=1, pch=pch,lwd=1.2)
    axis(2, at= seq(0,0.7,by=0.1), las=1,cex.axis = 1.0)
    
    dev.off()
    
    #compare = cbind(, rep(1, 4))
    #compare = rbind(compare, cbind(, rep(2, 2)))
    #compare = rbind(compare, cbind(NA, rep(3, 2)))
    #compare = rbind(compare, cbind(NA, rep(4, 2)))
    #compare = rbind(compare, cbind(, rep(5, 4)))
    #compare = rbind(compare, cbind(, rep(6, 2)))
    #mean(as.numeric(compare[which(compare[,2]==1), 1]))
    #boxplot(compare[,1]~compare[,2], ylim=c(0, 0.2))
    
    #### compare cell sizes
    x1 = keep$cell.size.mono.mean[c(4,12,20,28)]
    x11 = keep.ko$cell.size.mono.mean[c(1, 3)]
    x2 = keep$cell.size.binu.mean[c(4,12,20,28)]
    x22 = keep.ko$cell.size.binu.mean[c(1, 3)]
    x3 = keep$cell.size.mono.mean[c(8,16,24,32)]
    x33 = keep.ko$cell.size.mono.mean[c(2, 4)]
    x4 = keep$cell.size.binu.mean[c(8,16,24,32)]
    x44 = keep.ko$cell.size.binu.mean[c(2, 4)]
    
    t.test(x1, x11)
    t.test(x2, x22)
    t.test(x3, x33)
    t.test(x4, x44)
    t.test(x1, x3)
    
    ym = c(mean(keep$cell.size.mono.mean[c(4,12,20,28)]), mean(keep.ko$cell.size.mono.mean[c(1, 3)]),
            mean(keep$cell.size.binu.mean[c(4,12,20,28)]), mean(keep.ko$cell.size.binu.mean[c(1, 3)]), NA,
            mean(keep$cell.size.mono.mean[c(8,16,24,32)]), mean(keep.ko$cell.size.mono.mean[c(2, 4)]),
            mean(keep$cell.size.binu.mean[c(8,16,24,32)]), mean(keep.ko$cell.size.binu.mean[c(2, 4)])
            )
            
    dy = c(sd(keep$cell.size.mono.mean[c(4,12,20,28)])/sqrt(4), sd(keep.ko$cell.size.mono.mean[c(1, 3)])/sqrt(2),
            sd(keep$cell.size.binu.mean[c(4,12,20,28)])/sqrt(4), sd(keep.ko$cell.size.binu.mean[c(1, 3)])/sqrt(2), NA,
            sd(keep$cell.size.mono.mean[c(8,16,24,32)])/sqrt(4), sd(keep.ko$cell.size.mono.mean[c(2, 4)])/sqrt(2),
            sd(keep$cell.size.binu.mean[c(8,16,24,32)])/sqrt(4), sd(keep.ko$cell.size.binu.mean[c(2, 4)])/sqrt(2)
            )
    cols = c('gray70', 'gray30', 'gray70', 'gray30','white', 'gray70', 'gray30','gray70', 'gray30')
    
    pdf('myplots/FIGURES/cell_size_WT_KO.pdf', width=2.0, height=1.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(2,3,0.5,0.8)+0.1, tcl = -0.3)
    
    lims = range(c((ym-dy), (ym+dy)), na.rm=TRUE)
    barx = barplot(ym, col=cols, ylim=c(0, lims[2]), axes=FALSE, w=0.2, lwd=0.01)
    arrows(barx, ym-dy, barx, ym+dy, length=0.04, angle=90, code=3, col='black', lty=1, pch=pch,lwd=1.2)
    axis(2, at= c(0, 2000, 4000), las=1,cex.axis = 1.0)
    
    dev.off()
    
    ########
    #### Distribution of cell size and nuclear size and parameter estimation
    ########
    Plot.cell.nuclear.size = FALSE
    if(Plot.cell.nuclear.size)
    {
        source('functions_images.R')
        nucleus.all.fitting = TRUE
        cell.fitting = FALSE
        calibration.facs = TRUE
        Method = 'Normal'
        Filter.cells = TRUE
        nucleus.size.mean = FALSE
        scale = 0.34*0.34
        
        WT = TRUE
        
        if(WT)
        {
            dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
            nn = c(4, 8)
            NN = c(1:length(nn));
            
            plot.version = '_WT'
            
        }else{
            dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
            dna = dna[, c(1, 2, 4, 3, 5)]
            
            nn = c(9, 33, 21, 45)
            NN = c(1:length(nn))
            
            plot.version = '_KO'
        }
        
        for(n in NN)
        {
            
            if(WT) {
                tt = (nn[n]-1)*3;
            }else{
                tt = nn[n];
            }
            
            cat('ZT ', tt, '\n');
            test = c()
            
            if(WT)
            {
                m = 1;
                if(nn[n]<10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_0', nn[n], '_0', m, '.txt', sep='')
                if(nn[n]>=10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_', nn[n], '_0', m, '.txt', sep='')
            }else{
                filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/CRYDKO_try_2nd/Image_', nn[n], '_01', '_KO.txt', sep='')
            }
            
            aa = read.table(filename, sep='\t', header=FALSE)
            test = c(test, aa[1,2], aa[1, 3])  ### surface of image
            
            aa = aa[-1, ]
            colnames(aa) = c('cell.index', 'cell.size', 'nucleus.cell')
            
            ### remove cells with 3 nuclei
            bb = cbind(aa, find.repeating(aa[,1]))
            bb = bb[which(bb[,4]<3),]
            
            #### filter outliers of cells
            if(Filter.cells)
            {
                source('functions_images.R')
                cat('filter outlier cells \n')
                cc = cell.nucleus.filtering(bb)
                cc = data.frame(cc)
                colnames(cc) = c('cell.index', 'cell.size', 'nuclear.size.s', 'nb.nuclei', 'outlier')
                
                index.outlier = cc[which(cc[,5]==1) ,1]
                mm = match(bb[,1], index.outlier)
                mm = which(is.na(mm)==TRUE)
                bb = bb[mm, ]
                cc = cc[-which(cc[,5]==1), ]
            }
            
            #nucleus.size.mean = TRUE
            if(nucleus.size.mean){bb = cc;}
            
            cat(length(unique(bb[,1])), ' cells detected\n')
            
            test = c(test, length(unique(bb[,1])))
            
            percent = length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1)))
            test = c(test, length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1))))
            
            ### distribution of cell size
            jj = which(bb[,4]==1)
            kk = which(bb[,4]==2)
            
            cc = bb[,2]*scale
            xlim = c(0, 800)
            #ylim = c(0, 0.007)
            breaks = seq(0, 800, by=20)
            pdf.name = paste("myplots/FIGURES/Distribution_Mono_cell_size_ZT", tt, plot.version, ".pdf", sep='')
            pdf(pdf.name, width=1.7, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(0.25,0.25,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            hist(cc[jj], breaks=breaks, xlim=xlim, col='blue',freq=TRUE, xlab=NA, main=NA, ylab=NA, axes=FALSE, lwd=1.0)
            abline(v=median(cc[jj]), col='gray', lwd=1.5)
            
            axis(1,at=seq(0, 800, by=200),cex.axis =0.7)
            #axis(2, at= seq(0, 2000, by=400), las=0,cex.axis = 0.7)
            axis(2, at=c(0, 600, 1200), las=0,cex.axis = 0.7)
            box()
            
            dev.off()
            
            breaks = seq(0, (max(unique(cc[kk]))+20), by=20)
            pdf.name = paste("myplots/FIGURES/Distribution_Binu_cell_size_ZT", tt, plot.version, ".pdf", sep='')
            pdf(pdf.name, width=1.7, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(0.25,0.25,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            index = match(unique(bb[kk,1]), bb[kk, 1])
            
            hist(cc[kk[index]], breaks=breaks, xlim=xlim, freq=TRUE, col='red', xlab=NA, ylab=NA, main=NA, lwd=1.0, axes=FALSE)
            #abline(v=median(unique(cc[kk])), col='gray', lwd=1.5)
            abline(v=median(cc[kk[index]]), col='gray', lwd=1.5)
            
            axis(1,at=seq(0, 800, by=200),cex.axis =0.7)
            #axis(2, at= seq(0, 150, by=40), las=0,cex.axis = 0.7)
            axis(2, at=c(0, 50, 150, 250), las=0,cex.axis = 0.7)
            box()
            
            dev.off()
            
            #########################
            ### nuclear size analysis
            #########################
            ## all nucleui
            x = bb[,3];
            mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1))
            
            yy = x[which(x>250)]
            res = size.seperation(yy, method=Method, k=3, mu.init=mu.init)
            lambda = res[1:3]
            par1 = res[4:6]
            par2 = res[7:9]
            #ratio = c(par1[2]/par1[1], par1[3]/par1[2])
            refs = dna[c(1, 2, 4), (n+1)]
            alpha = refs/sum(refs)
            res = normalmixEM.test(x, k=3, mu.init=par1, sigma.init=par2, lambda.constr=alpha)
            par1 = res$mu
            
            if(n==31){
                cutoff.size=700
            }else{
                cutoff.size=600
            }
            nitr = 0
            while(par1[3]<cutoff.size)
            {
                mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
                res = normalmixEM.test(x, k=3, mu.init=mu.init, lambda.constr=alpha)
                par1 = res$mu
                cat(par1, '\n');
                nitr = nitr+1;
                if(nitr>=20) break;
            }
            
            res1 = c(res$lambda, res$mu, res$sigma)
            lambda = res1[1:3]
            par1 = res1[4:6]
            par2 = res1[7:9]
            par10 = par1
            par20 = par2
            
            ## nuclei within mononucleated cells
            x = bb[jj,3]
            res2 = size.seperation(x, method=Method, k=3, m.constr=par10, sigma.constr=par20)
            #res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
            lambda = res2[1:3]
            par1 = res2[4:6]
            par2 = res2[7:9]
            
            xfit<-seq(min(x),max(x),length=100)*scale
            yfit1<-dnorm(xfit, par1[1]*scale, par2[1]*scale)*lambda[1]
            yfit2<-dnorm(xfit, par1[2]*scale,par2[2]*scale)*lambda[2]
            yfit3<-dnorm(xfit, par1[3]*scale,par2[3]*scale)*lambda[3]
            
            xlim = c(10, 120)
            ylim = c(0, 0.04)
            #x = x*scale
            #xfit = xfit*scale
            pdf.name = paste("myplots/FIGURES/Distribution_Mono_Nuclear_size_ZT", tt, plot.version, ".pdf", sep='')
            pdf(pdf.name, width=1.7, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(0.25,0.25,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            breaks = seq(0, max(x*scale+10), by=3)
            hist(x*scale, breaks=breaks, xlim = xlim, ylim=ylim, xlab=NA, ylab=NA, freq=FALSE, col='blue', main=NA, axes=FALSE)
            lines(xfit, yfit1, col='green', lwd=1.2)
            lines(xfit, yfit2, col='green', lwd=1.2)
            lines(xfit, yfit3, col='green', lwd=1.2)
            #abline(v=par1*scale, col='green', lwd=1.2, lty=2)
            
            axis(1,at=c(10, 50, 100),las = 1, cex.axis =0.7)
            axis(2, at= seq(0, 0.04, by=0.02), las=0,cex.axis = 0.7)
            box()
            
            dev.off()

            ## nuclei within binucleated cells
            x = bb[kk,3]
            res3 = size.seperation(x, method=Method, k=3, m.constr=par10, sigma.constr=par20)
            
            lambda = res3[1:3]
            par1 = res3[4:6]
            par2 = res3[7:9]
            lambda.b = lambda
            
            xfit<-seq(min(x),max(x),length=100)*scale
            yfit1<-dnorm(xfit, par1[1]*scale, par2[1]*scale)*lambda[1]
            yfit2<-dnorm(xfit, par1[2]*scale, par2[2]*scale)*lambda[2]
            yfit3<-dnorm(xfit, par1[3]*scale,par2[3]*scale)*lambda[3]
            #test = c(test, percent*lambda[c(1:2)])
            
            pdf.name = paste("myplots/FIGURES/Distribution_Binu_Nuclear_size_ZT", tt, plot.version, ".pdf", sep='')
            pdf(pdf.name, width=1.7, height=1.6)
            par(cex = 0.7, las = 1, mgp = c(0.25,0.25,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
            
            breaks = seq(0, max(x*scale+10), by=3)
            hist(x*scale, breaks=50, xlim = xlim, ylim=ylim, xlab=NA, ylab=NA, freq=FALSE, col='red', main=NA, axes=FALSE)
            lines(xfit, yfit1, col='green', lwd=1.2)
            lines(xfit, yfit2, col='green', lwd=1.2)
            #abline(v=par1, col='green', lwd=1.2, lty=2)
            lines(xfit, yfit3, col='green', lwd=1.5)
            
            axis(1,at=c(10, 50, 100),las = 1, cex.axis =0.7)
            axis(2, at= seq(0, 0.04, by=0.02), las=0,cex.axis = 0.7)
            
            box()
            
            dev.off()
        
        }
    }
    
    ##################
    ##################
    #### To estimate parameters for S phase s(t), nuclear division n(t), mitosis m(t) by fitting percentages of cell populaitons
    ##################
    ##################
    Model.Cell.Cycle = FALSE
    if(Model.Cell.Cycle)
    {
        load(file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32.Rdata')
        remv = c(1, 26)
        kkeep = keep
        kkeep[remv, ] = NA
        
        fitting = FALSE
        if(fitting)
        {
            source('functions_images.R')
            
            percents = t(keep[-remv, c(14, 15, 17, 16, 18)])
            tt = setdiff(c(1:32), remv)*3
            
            res = fitting.polyploidy.cycle(time=tt, percents = percents)
        }
        
        load(file='Rdata/cell_endoreplication_fitting_model_selection.Rdata')
        
        o1 = order(rres$BIC)
        rres = rres[o1, ]

    }
    
}




