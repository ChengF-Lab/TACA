library(data.table);library(dplyr);library(ggplot2);library(ggrepel);library(mgcv)
library(ggh4x,lib='/home/lorincn/Rpkgs')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PLCG2 Microglia gene-specific Age-expression density plots ####
# data=readRDS('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/PLCG2_microglia_plot.Rds')
data=readRDS('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/PLCG2_microglia_plot.Rds')
plotdf=data$plotdf
plotdf$new_SampleID=rownames(plotdf)
meandf=data$meandf
corrdf=data$corrdf
gene.='PLCG2'
plotdf %>%
  filter(group %in% c('AD','PART','Resilience')) %>%
  filter(celltype %in% c('Inflammation','Tau','TCell')) %>%
  ggplot(aes(x=age,y=log2(1+y))) +
  # geom_point() +
  stat_density2d_filled(n=100,
                        # alpha=5/10,
                        bins=5,
                        contour_var="ndensity") +
  geom_hline(aes(yintercept=ybar),color='black',lwd=1/5,linetype='dotted',#data=meandf) +
             data=meandf %>%
               filter(group %in% c('AD','PART','Resilience')) %>%
               filter(celltype %in% c('Inflammation','Tau','TCell'))) +
  scale_fill_manual(values=colorRampPalette(c('white','#E7C3FA','#9C00F0'))(5)) +
  facet_grid(celltype~group) +
  stat_smooth(method='lm',color='black',fill='gray80',lwd=1/3,fullrange=TRUE) +
  labs(title=NULL,x='Age',y='Normalized gene expression of PLCG2') +
  theme_bw() +
  coord_cartesian(ylim=c(0,5)) +
  # lims(y=c(0,2)) +
  theme(legend.position='none',
        text=element_text(size=9,family="Arial"),
        axis.title=element_text(size=9,family='Arial'),
        axis.text=element_text(size=8,family='Arial'),
        strip.background=element_rect(fill="gray95"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=c(60,80,100),labels=c(60,80,100)) +
  geom_text(aes(x=85,y=4,label=label),size=2.35,#data=corrdf)
            data=corrdf %>%
              filter(group %in% c('AD','PART','Resilience')) %>%
              filter(celltype %in% c('Inflammation','Tau','TCell')))
setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave(paste0(gene.,'_ADvsHC_microglia_subtypes.svg'),width=3.25,height=3.25,units='in')
### testing mean differences across celltypes in AD vs Normal
B=plotdf %>% filter(group %in% c('AD','Normal')) %>% group_by(celltype,group) %>% summarise(xbar=mean(log2(1+y)),sbar=sd(log2(1+y))/sqrt(n())) %>% ungroup() %>% select(-c(sbar)) %>% tidyr::pivot_wider(names_from=group,values_from=xbar) %>% select(-celltype) %>% as.matrix()
S=plotdf %>% filter(group %in% c('AD','Normal')) %>% group_by(celltype,group) %>% summarise(xbar=mean(log2(1+y)),sbar=sd(log2(1+y))/sqrt(n())) %>% ungroup() %>% select(-c(xbar)) %>% tidyr::pivot_wider(names_from=group,values_from=sbar) %>% select(-celltype) %>% as.matrix()
Ds=(B[,1]-B[,2])/sqrt(rowSums(S^2))
stat=mean(Ds)^2/sqrt(sum(rowSums(S^2)))
D=(B[,1]-B[,2])^2/rowSums(S^2)
pchisq(D,1,lower.tail=FALSE)

### testing in just select datasets/cohorts
rdf=readRDS('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/metadata/Immune.Rds')
rdf=rdf %>% select(new_SampleID,Dataset) %>% distinct()
plotdf=plotdf %>% left_join(rdf,by='new_SampleID')

celltypes=unique(plotdf$celltype)
datasets=c('ROSMAP','SEA-AD')
resdf=data.frame()
for(i in 1:length(celltypes)) {
  for(j in 1:length(datasets)) {
    datij=plotdf %>%
      filter(celltype==celltypes[i],grepl(datasets[j],Dataset,fixed=TRUE)) %>%
      select(age,y) %>%
      as.matrix() %>%
      na.omit()
    datij[,2]=log2(datij[,2]+1)
    if(nrow(datij)<10) next
    r=cor(datij,method='spearman',use='complete.obs')[1,2]
    pval=cor.test(datij[,1],datij[,2],method='spearman',use='complete.obs')$p.value
    toadd=data.frame(celltype=celltypes[i],Dataset=datasets[j],r,pval,n=nrow(datij))
    resdf=rbind(resdf,toadd)
  }
}
resdf
d=sum(sign(resdf$r)==-1)
pbinom(d,nrow(resdf),1/2,lower.tail=FALSE)
table(resdf$r<0)
hist(resdf$r,breaks=20)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Microglia subtype volcano plot ####
# setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype_and_Homeostasis')
setwd('/mnt/isilon/w_gmi//chengflab/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype_and_Homeostasis')
data=readRDS('DEGs_ttest.Rds') %>% group_by(New_Group,target_celltype) %>% mutate(fdr=p.adjust(pval,'BH')) %>% ungroup()

set.seed(10)
data %>%
  # filter(target_celltype %in% c('Macrophage','TCell')) %>%
  filter(target_celltype %in% c('Tau','TCell')) %>%
  mutate(pval=ifelse(pval==0,min(pval[pval>0]),pval)) %>%
  mutate(mean_diff=mean_target-mean_homeostasis) %>%
  group_by(target_celltype,New_Group) %>%
  mutate(fdr=p.adjust(pval,'BH'),
         fdrsig=fdr<0.05 & abs(mean_diff)>0) %>%
  ungroup() %>%
  mutate(fdrlevel=cut(fdr,breaks=c(1,0.05,0.001,0))) %>%
  ggplot(aes(x=mean_diff,y=-log10(pval),color=factor(fdrlevel))) +
  geom_point(alpha=1/10,size=1/2) +
  scale_color_manual(values=c('#A357FAE4','#FC8B62','gray80')) +
  facet_grid2(target_celltype~New_Group,scales='free',independent='all') +
  theme_bw() +
  theme(strip.background=element_rect(fill="gray95")) +
  theme(legend.text=element_text(size=11),
        legend.position='none',
        text=element_text(size=9,family="Arial"),
        axis.title=element_text(size=9,family='Arial'),
        axis.text=element_text(size=6,family='Arial'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='mean difference in gene expression from Homeostasis',
       y=expression('-log'[10]*'(FDR)'))

setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('microglia_subtypes_WithinGroupBetweenCell_vs_HomeostasisOnly_DEGs.svg',width=4.5,height=2.5,units='in')
ggsave('microglia_subtypes_WithinGroupBetweenCell_vs_HomeostasisOnly_DEGs.pdf',width=4.5,height=2.5,units='in')
ggsave('microglia_subtypes_WithinGroupBetweenCell_vs_HomeostasisOnly_DEGs.png',width=4.5,height=2.5,units='in',dpi=1000)
# ~~~~~~~~~~~~~~~ #
# HPC !!!
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
# verifying association exists in presence of APOE status
setwd('/home/lorincn/isilon/Cheng-Qiu/TACA2/Integration_test/output_TACA2/Data_shared_Noah11/Meta/Immune')
refdata=fread('metadata.tsv.gz') %>% as.data.frame()
refdata=refdata %>% select(new_SampleID,APOE) %>% distinct() %>% filter(APOE!='Unknown')
# linear regression of plcg2 on age controlling for APOE
groups=c('AD','Resilience','PART')
cts=c('Inflammation','Tau','TCell')
rdf=data.frame()
for(j in 1:length(cts)) {
  setwd('~/isilon/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data/Immune')
  fp=paste0('all_cohorts_',cts[j],'_subtype_37824655_37824663_droppedCohorts_R6415047_droppedSubjects.Rds')
  df=readRDS(fp)
  metadata=df$pb_metadata
  counts=df$pb_counts
  plcg2=counts[which(rownames(counts)=='ENSG00000197943'),]/metadata$ncells
  plcg2=log2(1+plcg2)
  fitdf=bind_cols(metadata,data.frame(y=plcg2)) %>% inner_join(refdata,by='new_SampleID')
  fitdf=fitdf %>% mutate(Age=as.numeric(Age)) %>% tidyr::drop_na(Age)
  fitdf=fitdf %>%
    group_by(New_Group) %>%
    mutate(Age=impageall(Age))
  for(i in 1:length(groups)) {
    fit=lm(y~Age+APOE,data=fitdf %>% filter(New_Group==groups[i]))
    est=coef(fit)[2]
    se=sqrt(vcov(fit)[2,2])
    pval=pchisq((est/se)^2,1,lower.tail=FALSE)
    toadd=data.frame(group=groups[i],celltype=cts[j],est=est,se=se,pval=pval)
    rdf=rbind(rdf,toadd)
  }
}



# ~~~~~~~~~~~~~~~ #
# plot ADSP genes in AD group
adsp=read.csv('/Volumes/chengflab/Cheng-Noah/reference_data/ADSP_AD_genes.csv')
data %>%
  filter(symbol %in% adsp$gene) %>%
  filter(New_Group=='AD') %>%
  filter(!(target_celltype %in% c('Macrophage'))) %>%
  mutate(pval=ifelse(pval==0,min(pval[pval>0]),pval)) %>%
  mutate(mean_diff=mean_target-mean_homeostasis) %>%
  group_by(target_celltype,New_Group) %>%
  mutate(fdr=p.adjust(pval,'BH')) %>%
  ungroup() %>%
  ggplot(aes(x=mean_diff,y=symbol,color=fdr<0.05)) +
  geom_point(size=3/5) +
  geom_vline(xintercept=0,lwd=1/3) +
  facet_wrap(~target_celltype,nrow=1) +
  theme_bw() +
  scale_color_manual('',labels=c('FDR>0.05','FDR<0.05'),values=c('gray60','red')) +
  theme(legend.position='bottom') +
  theme(strip.background=element_rect(fill="gray95")) +
  labs(x='mean difference in gene expression from Homeostasis',
       y=NULL) +
  theme(axis.text.y=element_text(size=5),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.text=element_text(size=10)) +

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Sst Neuron volcano plot ####
library(dplyr);library(ggplot2);library(ggrepel)
# data=readRDS('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/Sst_Neuron_volcano.Rds')
data=readRDS('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/Sst_Neuron_volcano.Rds')
data=data %>%
  select(-fdr) %>%
  group_by(type,comparison) %>%
  mutate(fdr=p.adjust(p,'BH')) %>%
  ungroup()
data %>%
  filter(type=='tlm_withage') %>%
  filter(comparison=='AD vs Resilience') %>%
  filter(est<0,fdr<0.05) %>%
  filter(!grepl('MT-',symbol,fixed=TRUE)) %>%
  pull(est) %>%
  mean()
mt_genes=c('MT-ATP6',paste0('MT-CO',2:3),paste0('MT-ND',1:4),'MT-CYB')
data %>%
  filter(type=='tlm_withage') %>% # model type
  filter(comparison %in% c('AD vs Resilience')) %>%
  mutate(comparison=factor(comparison),
         comparison=recode_factor(comparison,
                                  'AD vs Resilience'='Resilience',
                                  'AD vs PART'='PART',
                                  'AD vs Normal'='Normal')) %>%
  mutate(fdrlevel=cut(fdr,breaks=c(1,0.05,0.001,0))) %>%
  ggplot(aes(x=est,y=-log10(fdr),color=fdrlevel)) +
  geom_point(alpha=1/4) +
  scale_color_manual(values=c('#A357FAE4','#FC8B62','gray80')) +
  theme_bw() +
  theme(strip.background=element_rect(fill="gray95")) +
  theme(text=element_text(size=8,family="Arial"),
        axis.title=element_text(size=8,family='Arial'),
        axis.text=element_text(size=7,family='Arial'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position='none') +
  scale_x_continuous(breaks=c(-0.8,0,0.7),labels=c(-0.7,0,0.7),limits=c(-1.,0.7)) +
  labs(title='AD vs Resilience',
       x=expression('-log'[2]*' effect size'),
       y=expression('-log'[10]*'(FDR q-value)')) +
  geom_point(aes(x=est,y=-log10(fdr)),color='firebrick',pch=19,size=1,data=. %>% filter(symbol %in% mt_genes)) +
  geom_text_repel(aes(x=est,y=-log10(fdr),label=symbol),
                  max.overlaps=100,
                  # nudge_x=-0.04,
                  force=2,
                  size=2,
                  segment.size=0.1,
                  # min.segment.length=0,
                  box.padding=0.1,
                  color='firebrick',
                  data=. %>% filter(symbol %in% mt_genes))
# setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('Sst_Neuron_DEGs_AD_Resilience_PART_Normal_tlm_withage.svg',width=1.9,height=1.4,units='in')
ggsave('Sst_Neuron_DEGs_AD_Resilience_PART_Normal_tlm_withage.pdf',width=1.9,height=1.4,units='in')
ggsave('Sst_Neuron_DEGs_AD_Resilience_PART_Normal_tlm_withage.png',width=1.9,height=1.4,units='in',dpi=700)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# MT-CYB western blot data from Tian ####
library(dplyr);library(ggplot2)
# data=read.csv('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/MT-CYB_western_blot.csv')
data=read.csv('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots/data/MT-CYB_western_blot.csv')
# corrdf=data.frame(path=c('ab','tau'),r=c('0.35','0.60'),p=c('5E-2','3E-3'))
data %>%
  filter(gene=='cyb') %>%
  ggplot(aes(gene_normalized,path_normalized)) +
  facet_wrap(~path) +
  stat_smooth(method='lm',se=FALSE,lwd=1,fill='gray90',color='firebrick') +
  coord_cartesian(ylim=c(0,220)) +
  geom_point(pch=19,color='#AF7AF0') +
  # geom_text(aes(x=60,y=200,label=paste0('r=',r,' (P=',p,')')),data=corrdf) +
  theme_bw() +
  theme(strip.background=element_rect(fill="gray95"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='CYB expression',y='protein abundance')
# setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
setwd('/mnt/isilon/w_gmi/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('CYB_western_blot_scatterplot_results.svg',width=3.75,height=2.1/5*4,units='in')
# correlations
library(magrittr)
bootci=function(data,gene.,path.,k=1000) {
  dat=data %>% filter(gene==gene.,path==path.)
  x=dat$gene_normalized
  y=dat$path_normalized
  n=nrow(dat)
  rs=c()
  for(i in 1:k) {
    ix=sample(1:n,n,replace=T)
    rs[i]=cor(x[ix],y[ix])
  }
  rs
}
abeta_cor=bootci(data,'cyb','amyloid-beta',1e4)
tau_cor=bootci(data,'cyb','pTau231',1e4)
mean(abeta_cor)
mean(tau_cor)
1-mean(abeta_cor>0)
1-mean(tau_cor>0)
#

data %>%
  filter(gene=='cyb') %>%
  ggplot(aes(y=gene_normalized,)) +
  geom_boxplot() +
  facet_wrap(~path)


#


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Global Age-Expression correlations in AD and Normal groups for microglia subtypes ####
library(dplyr);library(tidyr);library(ggplot2)
kernel=function(x,y,H=NULL) {
  X=cbind(x,y)
  if(is.null(H)) H=ks::Hpi(X)
  E=eigen(H)
  E=E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
  n=nrow(X)
  Z=matrix(nr=n,nc=n)
  h=sqrt(det(H))
  k=c()
  for(i in 1:n) {
    newXY=matrix(X[i,],nr=nrow(X),nc=2,byrow=T)
    k[i]=1/(nrow(X)*h)*sum(dmvn((newXY-X)%*%E,c(0,0),diag(2)))
  }
  k
}
gridkernel=function(grid1,grid2,x,y,H=NULL) {
  X=cbind(x,y)
  if(is.null(H)) H=ks::Hpi(X)
  E=eigen(H)
  E=E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
  GRID=cbind(grid1,grid2)
  ngrid=length(grid1)
  nobs=length(x)
  adj=1/(n*det(H))
  k=c()
  for(i in 1:ngrid) {
    Gi=matrix(GRID[i,],nrow(X),2,byrow=T)
    diff=Gi-X
    d=mvnfast::dmvn(diff%*%E,c(0,0),diag(2))
    k[i]=sum(d)
  }
  k
}
plotdf=readRDS('/Users/lorincn/Documents/microglia_subtype_plotdf.Rds')
plotdf=plotdf %>% filter(fdr<0.05)
cts=unique(plotdf$celltype)
pdf=regiondf=textdf=data.frame()
for(ct in 1:length(cts)) {
  XY=plotdf %>% 
    filter(celltype==cts[ct]) %>%
    select(AD,Normal) %>%
    as.matrix()
  H=ks::Hpi(XY)
  n=nrow(XY)
  x=XY[,1]
  y=XY[,2]
  gridn=1000
  newx=seq(min(x)*0.9,max(x)/0.9,length.out=gridn)
  newy=seq(min(y)*0.9,max(y)/0.9,length.out=gridn)
  grid=expand.grid(newx,newy)
  kd=gridkernel(grid[,1],grid[,2],x,y,H)
  kdf=data.frame(x=grid$Var1,y=grid$Var2,kd=kd)
  Z=matrix(kd,nr=gridn,nc=gridn)
  # plot_ly(x=newx,y=newy,z=Z,colors=c('#EBCC2A','#F21A00')) %>% add_surface()
  stdkd=kd/sum(kd)
  kdf$stdkd=stdkd
  M0ix=which.max(stdkd)
  M0=c(grid[M0ix,1],grid[M0ix,2])
  eval_zs=unique(round(stdkd,5))
  vol=c()
  for(i in 1:length(eval_zs)) {
    kix=which(round(kdf$stdkd,5)>=eval_zs[i])
    cubes=kdf[kix,] %>% select(x,y,stdkd)
    vol[i]=sum(cubes$stdkd)
  }
  c0=0.5
  wv=which.min(abs(vol-c0))
  z0=eval_zs[wv]
  stdZ=Z/sum(c(Z))
  GRIDcut=as.matrix(grid)[which(round(stdkd,5)>=z0),]
  # save
  pdf=rbind(pdf,data.frame(x=XY[,1],y=XY[,2],celltype=cts[ct]))
  hull_indices=chull(GRIDcut)
  hull_indices=c(hull_indices,hull_indices[1]) # close loop if its open
  Gch=GRIDcut[hull_indices,]
  regiondf=rbind(regiondf,data.frame(x=Gch[,1],y=Gch[,2],celltype=cts[ct]))
  textdf=rbind(textdf,data.frame(group=c('AD','Normal'),r=M0,celltype=cts[ct]))
}
textdf=textdf %>% pivot_wider(values_from='r',names_from='group')
textdf=textdf %>%
  mutate(AD_label=paste0(round(AD,2)),
         Normal_label=paste0(round(Normal,2)))
pdf %>%
  filter(celltype %in% c('DAM','MHC','TCell','Macrophage')) %>%
  ggplot(aes(x,y)) +
  facet_wrap(~celltype) +
  geom_point(col='white') +
  geom_hline(yintercept=0,color='gray80',lwd=1/3) +
  geom_vline(xintercept=0,color='gray80',lwd=1/3) +
  geom_segment(aes(x=-1,xend=AD,y=Normal,yend=Normal),color='#8AEB94',lwd=1/3,
               data=textdf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  geom_segment(aes(x=AD,xend=AD,y=-1,yend=Normal),color='#8AEB94',lwd=1/3,
               data=textdf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  geom_text(aes(x=-0.425,y=Normal+0.05,label=Normal_label),size=2.5,color='#0C8518',fontface='bold',
            data=textdf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  geom_text(aes(x=AD-0.05,y=-0.4,label=AD_label),size=2.5,color='#0C8518',fontface='bold',angle=90,
            data=textdf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  geom_path(aes(x,y),color='#0C8518',
            data=regiondf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  geom_point(aes(AD,Normal),pch=19,color='#0C8518',size=1/3,
             data=textdf %>% filter(celltype %in% c('DAM','MHC','TCell','Macrophage'))) +
  coord_cartesian(xlim=c(-0.5,0.1),ylim=c(-0.6,0.25)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  theme(strip.background=element_rect(fill="gray95")) +
  labs(x='correlation in AD subjects',
       y='correlation in cognitively normal subjects')
setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('age_expression_correlation_AD_vs_Normal_microglia_subtypes.svg',width=3,height=3,units='in')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# trying to make a plot to illustrate immune aging clock (barplot) ####
library(dplyr);library(tidyr);library(ggplot2)
kernel=function(x,y,H=NULL) {
  X=cbind(x,y)
  if(is.null(H)) H=ks::Hpi(X)
  E=eigen(H)
  E=E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
  n=nrow(X)
  Z=matrix(nr=n,nc=n)
  h=sqrt(det(H))
  k=c()
  for(i in 1:n) {
    newXY=matrix(X[i,],nr=nrow(X),nc=2,byrow=T)
    k[i]=1/(nrow(X)*h)*sum(dmvn((newXY-X)%*%E,c(0,0),diag(2)))
  }
  k
}
gridkernel=function(grid1,grid2,x,y,H=NULL) {
  X=cbind(x,y)
  if(is.null(H)) H=ks::Hpi(X)
  E=eigen(H)
  E=E$vectors%*%diag(1/sqrt(E$values))%*%t(E$vectors)
  GRID=cbind(grid1,grid2)
  ngrid=length(grid1)
  nobs=length(x)
  adj=1/(n*det(H))
  k=c()
  for(i in 1:ngrid) {
    Gi=matrix(GRID[i,],nrow(X),2,byrow=T)
    diff=Gi-X
    d=mvnfast::dmvn(diff%*%E,c(0,0),diag(2))
    k[i]=sum(d)
  }
  k
}
plotdf=readRDS('/Users/lorincn/Documents/microglia_subtype_plotdf.Rds')
longrg=plotdf %>%
  select(AD,Normal,Resilience,PART,symbol,celltype,fdr) %>%
  distinct() %>%
  pivot_longer(cols=c('AD','Normal','Resilience','PART')) %>%
  rename(group=name,rg=value)
ordy=longrg %>%
  filter(group=='AD') %>%
  group_by(symbol) %>%
  mutate(avgrg=mean(rg)) %>%
  ungroup() %>%
  arrange(avgrg) %>%
  pull(symbol) %>%
  unique()
longrg %>%
  filter(celltype %in% c('Inflammation','TCell','Tau','Homeostasis')) %>%
  mutate(group=factor(group,levels=c('Normal','AD','PART','Resilience'))) %>%
  mutate(celltype=factor(celltype,levels=c('Tau','Homeostasis','TCell','Inflammation'))) %>%
  mutate(symbol=factor(symbol,levels=ordy)) %>%
  ggplot(aes(x=rg,y=symbol,color=rg<0)) +
  scale_color_manual(values=c('#7415E8','#10B31E')) +
  geom_bar(stat='identity') +
  facet_grid(celltype~group) +
  theme_bw() +
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1),labels=c(-1,'',0,'',1),limits=c(-1,1)) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  theme(strip.background=element_rect(fill="gray95")) +
  theme(legend.position='none',
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x='correlation between age and gene expression',
       y='gene')
setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/final_plots')
ggsave('age_expression_correlation_AD_vs_Normal_microglia_subtypes_barplot.svg',width=3.25,height=6,units='in')










