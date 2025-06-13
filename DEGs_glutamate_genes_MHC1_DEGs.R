# rm(list=ls(all=TRUE))
library(data.table);library(dplyr);library(ggplot2);
library(beeswarm)
library(ggbeeswarm,lib='/home/lorincn/Rpkgs')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
###
ens=c('ENSG00000196277','ENSG00000168959','ENSG00000273079','ENSG00000183454','ENSG00000152578','ENSG00000120251')
genes=c('GRM7','GRM5','GRIN2B','GRIN2A','GRIA4','GRIA2')
###
# MHC1 data
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data/Immune/')
data=readRDS('all_cohorts_MHC1_subtype_37824655_37824663_droppedCohorts_R6415047_droppedSubjects_99pctQC.Rds')
metadata=data$pb_metadata %>% mutate(Age=ifelse(Age=='Unknown',NA,Age)) %>% mutate(Age=as.numeric(Age))
counts=data$pb_counts
NC=matrix(metadata$ncells,nr=nrow(counts),nc=ncol(counts))
counts=counts[ens,]
countdf=data.frame(y=c(counts),
                   gene=rep(genes,ncol(counts)),
                   ensembl_id=rep(ens,ncol(counts)),
                   new_SampleID=rep(metadata$new_SampleID,each=nrow(counts)))
countdf=left_join(countdf,metadata,by='new_SampleID')
meandf=countdf %>%
  group_by(gene,New_Group) %>%
  summarise(x=mean(log2(1+y/ncells)),
            se=sd(log2(1+y/ncells))/sqrt(n()))
countdf %>%
  # filter(New_Group %in% c('AD','Resilience')) %>%
  filter(New_Group %in% c('AD','Normal')) %>%
  # filter(New_Group %in% c('AD','Resilience','Normal','PART')) %>%
  ggplot(aes(gene,log2(1+y/ncells),fill=New_Group)) +
  scale_fill_manual('',values=c('#EAB06E','#019A92')) +
  # geom_violin(position=position_dodge(3/4)) +
  geom_boxplot(width=3/6,position=position_dodge(3/4),outlier.size=1/2) +
  theme_bw() +
  coord_cartesian(ylim=c(0,2.5)) +
  theme(legend.position='bottom',
        # panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x=NULL,y='std. gene expression')
X_AD=meandf %>% filter(New_Group=='AD')
X_NO=meandf %>% filter(New_Group=='Normal')
stats=(X_AD$x-X_NO$x)^2/(X_AD$se^2+X_NO$se^2)
data.frame(gene=X_AD$gene,pvalue=pchisq(stats,1,lower.tail=FALSE)) %>%
  mutate(fdr=p.adjust(pvalue,'BH'))
