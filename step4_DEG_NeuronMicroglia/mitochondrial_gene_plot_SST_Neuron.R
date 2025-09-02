# rm(list=ls(all=TRUE))
library(data.table);library(dplyr);library(ggplot2)
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
  if(length(nn)==0) return(ages)
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
metadata=metadata %>% group_by(Sex,New_Group) %>% mutate(Age=impageall(Age)) %>% ungroup()

fit_data=metadata %>% filter(New_Group %in% c('AD','Normal')) %>% mutate(AD=+(New_Group=='AD'))
fit=glm(AD~I(Sex=='Male')+Age,data=fit_data,family=binomial(link='logit'))
summary(fit)
ages=range(fit_data$Age,na.rm=TRUE)
male=c(1,0)
df=expand.grid(male=male,age=seq(ages[1],ages[2],length.out=30))
y=cbind(1,as.matrix(df))%*%coef(fit)
y=1/(1+exp(-y))
df$y=c(y)
ggplot(df,aes(age,y,color=factor(male))) +
  geom_path() +
  scale_color_manual('',labels=c('Females','Males'),values=c('red','blue')) +
  theme_bw() +
  theme(legend.position='bottom',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='Age',y='P(AD|Age,Sex)')

#


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
# plotting select mitochondrial genes
plotdf %>%
  filter(New_Group %in% c('Normal','Resilience','AD')) %>%
  mutate(New_Group=factor(New_Group,levels=c('Normal','Resilience','AD'))) %>%
  ggplot(aes(New_Group,log2(1+y),fill=New_Group)) +
  geom_boxplot() +
  facet_wrap(~gene,nrow=2) +
  scale_fill_manual(values=c('white','#AF7AF0','#0C8518')) +
  theme_bw() +
  lims(y=c(0,10)) +
  theme(legend.position='none',
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  theme(strip.background=element_rect(fill="gray95")) +
  labs(x='disease group',y=expression('log'[2]*'(1+pseudobulk expression)'))
# mean expression in nuclear mitochrondrial genes
metadata=data$pb_metadata %>% mutate(Age=ifelse(Age=='Unknown',NA,Age)) %>% mutate(Age=as.numeric(Age))
counts=data$pb_counts
NC=matrix(metadata$ncells,nr=nrow(counts),nc=ncol(counts),byrow=TRUE)
counts=counts/NC
mitocarta_in=mitocarta[mitocarta %in% rownames(counts)]
counts=counts[mitocarta_in,]
ADix=which(metadata$New_Group=='AD')
REix=which(metadata$New_Group=='Resilience')
NOix=which(metadata$New_Group=='Normal')
yAD=rowMeans(log2(counts[,ADix]+1))
yRE=rowMeans(log2(counts[,REix]+1))
yNO=rowMeans(log2(counts[,NOix]+1))
# test all one by one
XAD=log2(counts[,ADix]+1)
XRE=log2(counts[,REix]+1)
D=rowMeans(XAD)-rowMeans(XRE)
# V=apply(XAD,1,mad)^2/ncol(XAD)+apply(XRE,1,mad)^2/ncol(XRE)
V=apply(XAD,1,var)/ncol(XAD)+apply(XRE,1,var)/ncol(XRE)
stats=D^2/V
ps=pchisq(stats,1,lower.tail=FALSE)
fdr=p.adjust(ps,'BH')
plot(-log10(sort(fdr)))
table(ps<0.05)
table(fdr<0.05)
table(ps<(0.05)/length(ps))
# mean across people for each gene
xbars=c(mean(yAD),
        mean(yRE),
        mean(yNO))
sbars=c(sd(yAD)/sqrt(length(ADix)),
        sd(yRE)/sqrt(length(REix)),
        sd(yNO)/sqrt(length(NOix)))
# Resilience vs AD
xbars[2]-xbars[1]
stat=(xbars[2]-xbars[1])^2/(sbars[2]^2+sbars[1]^2)
pchisq(stat,1,lower.tail=FALSE)
# Resilience vs Normal
xbars[2]-xbars[3]
stat=(xbars[2]-xbars[3])^2/(sbars[2]^2+sbars[3]^2)
pchisq(stat,1,lower.tail=FALSE)
# binomial test Resilience vs AD
d=+(yRE>yAD)
pbinom(sum(d),length(yAD),prob=1/2,lower.tail=FALSE)
# binomial test Resilience vs Normal
d=+(yRE>yNO)
pbinom(sum(d),length(yAD),prob=1/2,lower.tail=FALSE)

# nicer plot
df=data.frame(
  diff=c(yRE-yAD,yRE-yNO),
  type=rep(c('Resilience - AD','Resilience - Normal'),each=length(yRE)),
  gene=rep(rownames(counts),2)
)
df=df %>%
  group_by(type) %>%
  arrange(diff) %>%
  mutate(ind=1:n()) %>%
  ungroup()
(p1=df %>%
  filter(type=='Resilience - AD') %>%
  ggplot(aes(ind,diff,color=diff<0,alpha=abs(diff)>log2(1.1))) +
  geom_segment(aes(x=ind,xend=ind,y=0,yend=diff)) +
  scale_color_manual(values=c('#3C508B','#4FC46A')) +
  geom_hline(yintercept=c(1,-1)*log2(1.1),lwd=1/2,linetype='dashed') +
  geom_hline(yintercept=0,lwd=1/2) +
  theme_bw() +
  scale_alpha_manual(values=c(1/3,1)) +
  theme(legend.position='none',
        strip.background=element_rect(fill="gray95"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
    scale_y_continuous(breaks=seq(-0.4,0.4,0.1),labels=seq(-0.4,0.4,0.1)) +
  labs(x='sorted mitochondrial gene (nDNA)',
       y=expression('log'[2]*' fold change'),
       title='CR vs AD'))
(p2=df %>%
  filter(type=='Resilience - Normal') %>%
  ggplot(aes(ind,diff,color=diff<0)) +
  geom_segment(aes(x=ind,xend=ind,y=0,yend=diff)) +
  scale_color_manual(values=c('#3C508B','#4FC46A')) +
  theme_bw() +
  theme(legend.position='none',
        strip.background=element_rect(fill="gray95"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
    lims(y=c(-0.2,0.6)) +
  labs(x='mitochondrial gene index (DNA)',
       y='expression (CR) - expression(CN)',
       title='CR vs CN'))
ggpubr::ggarrange(p1,p2,nrow=1)
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('nuclear_mitochondrial_DNA_genomewide_plot.pdf',width=5.54,height=2.94)






































