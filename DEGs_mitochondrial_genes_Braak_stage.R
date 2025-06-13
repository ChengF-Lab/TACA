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
# nuclear mitochondria genes
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2')
mitocarta=read.csv('MitoCarta3_mitochrondrial_genes.csv') %>% 
  select(chr=hg19_Chromosome,ensembl_id=EnsemblGeneID_mapping_version_20200130) %>%
  filter(!(chr %in% c('','chrM','chrX'))) %>%
  pull(ensembl_id) %>%
  sapply(.,function(h) unlist(strsplit(h,'[|]')))
mitocarta=unlist(unname(mitocarta))
# Sst data
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data/Sst_Neuron')
data=readRDS('all_cohorts_Sst_Neuron_37824655_37824663_droppedCohorts_R6415047_droppedSubjects_99pct_Qc.Rds')
metadata=data$pb_metadata %>% mutate(Age=ifelse(Age=='Unknown',NA,Age)) %>% mutate(Age=as.numeric(Age))
counts=data$pb_counts
mt_genes=matrix(c(
  'MT-ATP6','ENSG00000198899',
  'MT-CO2','ENSG00000198712',
  'MT-CO3','ENSG00000198938',
  'MT-ND1','ENSG00000198888',
  'MT-ND2','ENSG00000198763',
  'MT-ND3','ENSG00000198840',
  'MT-ND4','ENSG00000198886',
  'MT-CYB','ENSG00000198727'),
  byrow=T,nc=2
)
colnames(mt_genes)=c('symbol','ensembl_id')
mt_genes=as.data.frame(mt_genes)
ix=which(rownames(counts) %in% mt_genes$ensembl_id)
counts=counts[ix,]
NC=matrix(metadata$ncells,nr=nrow(counts),nc=ncol(counts),byrow=TRUE)
counts=counts/NC
counts=counts[mt_genes$ensembl_id,]
countdf=data.frame(
  y=c(counts),
  new_SampleID=rep(metadata$new_SampleID,each=nrow(counts)),
  gene=rep(mt_genes$symbol,ncol(counts))
)
plotdf=left_join(countdf,metadata,by='new_SampleID')
# recode Braak Stage
plotdf$BraakStage=as.character(plotdf$BraakStage)
missing=c('Unknown','N')
zeros=c('0','0.5','0-I','0.0')
ones=c('1.0','1.5','I')
twos=c('2.0','2.5','II','II-III')
threes=c('3.0','3.5','III')
fours=c('4.0','IV')
fives=c('5.0','V')
sixes=c('VI')
plotdf$BraakStage[plotdf$BraakStage %in% missing]=NA
plotdf$BraakStage[plotdf$BraakStage %in% zeros]='0'
plotdf$BraakStage[plotdf$BraakStage %in% ones]='1'
plotdf$BraakStage[plotdf$BraakStage %in% twos]='2'
plotdf$BraakStage[plotdf$BraakStage %in% threes]='3'
plotdf$BraakStage[plotdf$BraakStage %in% fours]='4'
plotdf$BraakStage[plotdf$BraakStage %in% fives]='5'
plotdf$BraakStage[plotdf$BraakStage %in% sixes]='6'
plotdf$BraakStage=as.numeric(plotdf$BraakStage)
# plot
(cordf=plotdf %>%
    filter(BraakStage>0) %>%
  na.omit() %>%
  group_by(gene) %>%
  summarise(r=cor(log2(1+y),BraakStage),
            pval=cor.test(log2(1+y),BraakStage)$p.value))
plotdf %>%
  na.omit() %>%
  filter(gene %in% cordf$gene[cordf$pval<0.05]) %>%
  filter(gene %in% c('MT-CYB','MT-ATP6','MT-CO2')) %>%
  filter(BraakStage>0) %>%
  group_by(gene) %>%
  mutate(y=y/max(y)) %>%
  ungroup() %>%
  ggplot(aes(BraakStage,log2(1+y))) +
  geom_beeswarm(size=1/3,color='#7B00FF1A') +
  facet_wrap(~gene,nrow=1) +
  stat_smooth(method='lm',se=FALSE,color='#03A821') +
  theme_bw() +
  theme(strip.background=element_rect(fill="gray95"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='Braak stage',y='std. expression') +
  geom_text(aes(x=3.5,y=0.9,label=lab),
            size=2.85,
            data=cordf %>%
              filter(gene %in% c('MT-CYB','MT-ATP6','MT-CO2')) %>%
              filter(pval<0.05) %>%
              mutate(lab=paste0('r=',round(r,2),'\nP=',formatC(pval,format='e',digits=1))))
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('mitochondrial_gene_Braak_stage_correlation_maintext.svg',width=3,height=1.8)


















