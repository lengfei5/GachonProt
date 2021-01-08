##########################################
# localization annotation
##########################################
find.uniprot.localization = function(subs.u, annots.local = c('Nucleus', 'Endoplasmic reticulum', 'Golgi apparatus', 
                                                              'Cell membrane', 'Vesicle/Vacuole',
                                                              'Mitochondrion', 'Extracellular space', 'Cytosol', 
                                                              'Cytoskeleton', 'Peroxisome', 'Lysosome',  'Endosome'))
{
  find = c()
  for(p in 1:12)
  {
    #p = 1;
    kk = c(grep(annots.local[p], subs.u), grep(tolower(annots.local[p]), subs.u))
    if(p==1)  kk = c(kk, grep('Nucleolus', subs.u), grep(tolower('Nucleolus'), subs.u))
    if(p==4)  kk = c(kk, grep('Plasma membrane', subs.u), grep(tolower('Plasma membrane'), subs.u))
    if(p==5)  kk = c(kk, grep('Vesicle', subs.u), grep(tolower('Vesicle'), subs.u), grep('Vacuole', subs.u), grep(tolower('Vacuole'), subs.u))
    if(p==7)  kk = c(kk, grep('Extracellular region', subs.u), grep(tolower('Extracellular region'), subs.u))
    if(p==8)  kk = c(kk, grep('Cytoplasm', subs.u), grep(tolower('Cytoplasm'), subs.u), grep('Cytoplasmic part', subs.u), grep('cytoplasmic part', subs.u))
    
    if(length(kk)>0) find = c(find, p);
  }
  if(length(setdiff(c(8, 9), find))==0) find = setdiff(find, c(8))
  if(length(find)>0) 
  {
    find = find[order(find)];
    return(paste(annots.local[find], sep='', collapse = ';'));
  }else{
    return(NA);
  }
  
}
