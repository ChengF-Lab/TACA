library(dplyr);library(magrittr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions
le=function(x,y,x0) {
  # linearly extrapolate value between two points in 2D
  XY=cbind(x,y)
  XY=XY[order(XY[,1]),]
  ix=which(x>x0)[1]
  if(length(ix)==0) return(tail(XY[,2],1))
  xa=XY[ix-1,1]
  xb=XY[ix,1]
  ya=XY[ix-1,2]
  yb=XY[ix,2]
  m=(yb-ya)/(xb-xa)
  a=yb-xb*m
  y0=a+x0*m
  y0
}
# no need to impute
rkd=function(n,fullx,a,b=NA) {
  if(is.na(b)) b=max(fullx)
  d=density(fullx)
  x=d$x;y=d$y
  data=c()
  while(length(data)<n) {
    r1=runif(1,a,b)
    r2=runif(1,0,max(y))
    ytrue=le(x,y,r1)
    if(r2<ytrue) data=c(data,r1)
  }
  data
}
impage=function(ages,val=80.001) {
  aa=ages
  aa[aa==val]=NA
  ix=is.na(aa)
  l=sum(ix)
  bb=na.omit(aa)
  ia=rkd(l,bb,round(val))
  aa[ix]=ia
  aa
}
impageall=function(ages) {
  aa=as.character(ages)
  nn=table(aa[which(grepl('.001',aa))]) # doesn't have any '.999's
  nn=as.numeric(names(nn))
  nn=sort(nn)
  ix=which(as.numeric(ages) %in% nn)
  aa[ix]=NA
  aa=as.numeric(aa)
  for(i in 1:length(nn)) {
    ixi=which(as.numeric(ages) %in% nn[i])
    impage=rkd(length(ixi),na.omit(aa),nn[i],ifelse(is.na(nn[i+1]),NA,nn[i+1]))
    ages[ixi]=as.character(impage)
  }
  as.numeric(ages)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data/Immune/')
data=readRDS('all_cohorts_MHC1_subtype_37824655_37824663_droppedCohorts_R6415047_droppedSubjects_99pctQC.Rds')
metadata=data$pb_metadata
counts=data$pb_counts
rm(data)
ix=which(metadata$New_Group %in% c('AD','Normal'))
metadata=metadata[ix,]
counts=counts[,ix]
# impute age
metadata=metadata %>%
  group_by(New_Group) %>%
  mutate(Age=impageall(Age))
Xdf=metadata %>% mutate(AD=ifelse(New_Group=='AD',1,0)) %>% select(AD,Age,Sex)
# normalized by cell counts
NC=matrix(metadata$ncells,nrow(counts),ncol(counts),byrow=TRUE)
counts=counts/NC
genes=unique(rownames(counts))
rdf=data.frame()
for(i in 1:length(genes)) {
    if(i%%round(nrow(counts)*0.01)==0) cat(round(i*100/nrow(counts)),'% complete\n',sep='')
    fitdf=bind_cols(data.frame(y=counts[i,]),Xdf) %>% na.omit()
    fit=lm(I(log2(1+y))~Age+Sex+AD,data=fitdf)
    est=tail(coef(fit),1)
    se=tail(sqrt(diag(vcov(fit))),1)
    pval=pchisq((est/se)^2,1,lower.tail=FALSE)
    rsq=summary(fit)$r.squared
    toadd=data.frame(
        gene=genes[i],
        log2FC=est,
        SE=se,
        Pvalue=pval,
        pct_expressed_AD=mean(fitdf$y[fitdf$AD==1]>0),
        pct_expressed_Normal=mean(fitdf$y[fitdf$AD==0]>0),
        mean_ncells_AD=mean(NC[i,which(fitdf$AD==1)]),
        mean_ncells_Normal=mean(NC[i,which(fitdf$AD==0)])
    )
    rdf=rbind(rdf,toadd)
}
rdf=rdf %>% mutate(fdr_all=p.adjust(Pvalue,'BH')) %>% as_tibble()
rdf=rdf %>%
    mutate(istenpct=pct_expressed_AD>0.1 & pct_expressed_Normal>0.1) %>%
    group_by(istenpct) %>%
    mutate(fdr_tenpct=p.adjust(Pvalue,'BH')) %>%
    ungroup() %>%
    mutate(fdr_tenpct=ifelse(istenpct,fdr_tenpct,NA)) %>%
    select(-istenpct)
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/results/new/MHC1')
write.csv(rdf,'DEG_test_results_99pctQC.csv',row.names=FALSE,quote=FALSE)
# just DEGs, clean format
rdf %>% 
    filter(fdr_all<0.05) %>%
    select(gene,log2FC,fdr=fdr_all) %>%
    write.csv('DEGs_99pctQC.csv',row.names=FALSE,quote=FALSE)
rdf %>% 
    filter(fdr_tenpct<0.05) %>%
    select(gene,log2FC,fdr=fdr_all) %>%
    write.csv('DEGs_90pctQC.csv',row.names=FALSE,quote=FALSE)
