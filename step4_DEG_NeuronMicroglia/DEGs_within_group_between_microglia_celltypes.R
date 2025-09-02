# HPC!
library(data.table);library(dplyr);library(magrittr);library(tidyr);library(wesanderson)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions ####
cleanres=function(fit) {
  est=tail(coef(fit),1)
  se=tail(sqrt(diag(vcov(fit))),1)
  p=pchisq((est/se)^2,1,lower.tail=FALSE)
  rsq=summary(fit)$r.squared
  list(est=est,se=se,p=p,rsq=rsq)
}
cleanresdf=function(fit) as_tibble(cleanres(fit))
rtruncnorm=function(n,a,b,mu,sigma2) {
  x=c()
  i=0
  while(length(x)<n) {
    i=i+1
    xi=rnorm(n,mu,sqrt(sigma2))
    wb=which(xi>a & xi<b)
    if(length(wb)>0) x[i]=sample(xi[wb],1)
  }
  x
} # rejection sampling
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# analysis ####
groups=c('AD','Normal','PART','Resilience')
celltypes=c('DAM','Homeostasis','Inflammation','Macrophage','MHC','Tau','Proliferation','TCell')
# celltypes=celltypes[celltypes!='Homeostasis']
rdf=data.frame()
for(i in 1:length(groups)) {
  if(i==1) cat('\n')
  cat(groups[i],'\n')
  # load all data
  all_metadata=data.frame()
  all_counts=list()
  for(j in 1:length(celltypes)) {
    # load pseudobulk data
    setwd('~/isilon/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data/Immune')
    fpin=paste0('all_cohorts_',celltypes[j],'_subtype_37824655_37824663_droppedCohorts_R6415047_droppedSubjects.Rds')
    pb_dati=readRDS(fpin)
    metadata=pb_dati$pb_metadata
    counts=pb_dati$pb_counts
    ix=which(metadata$New_Group==groups[i])
    metadata=metadata[ix,]
    counts=counts[,ix]
    # ~~~ normalize by cell counts
    NC=matrix(metadata$ncells,nr=nrow(counts),nc=ncol(counts),byrow=TRUE)
    counts=counts/NC
    # ~~~ normalize by cell counts
    all_metadata=rbind(all_metadata,metadata %>% mutate(celltype=celltypes[j]))
    all_counts[[j]]=counts; names(all_counts)[j]=celltypes[j]
    if(j==1) genes=rownames(counts) else genes=intersect(genes,rownames(counts))
  }
  # put all data together into single objects
  all_counts_together=matrix(nr=length(genes),nc=1)
  celltype_indices=list()
  for(o in 1:length(all_counts)) {
    ci=all_counts[[o]]
    ix=which(rownames(ci) %in% genes)
    all_counts_together=cbind(all_counts_together,ci[ix,])
    if(o==1) all_counts_together=all_counts_together[,-1]
    celltype_indices[[o]]=(ncol(all_counts_together)-ncol(ci)+1):ncol(all_counts_together)
    names(celltype_indices)[o]=names(all_counts)[o]
  }
  # prepare data for single comparison (ie one celltype vs another)
  for(j in 1:length(celltypes)) {
    cat(' ',celltypes[j],'\n',sep='')
    ref=celltypes[j]
    nonref=celltypes[-j]
    counts_ref=all_counts_together[,celltype_indices[[j]]]
    metadata_ref=all_metadata %>% filter(celltype==celltypes[j])
    counts_nonref=all_counts_together[,-celltype_indices[[j]]]
    metadata_nonref=all_metadata %>% filter(celltype!=celltypes[j])
    # subset to just individuals with cells in all celltypes
    useppl=intersect(metadata_ref$new_SampleID,metadata_nonref$new_SampleID)
    metadata_ref=metadata_ref[which(colnames(counts_ref) %in% useppl),]
    counts_ref=counts_ref[,which(colnames(counts_ref) %in% useppl)]
    metadata_nonref=metadata_nonref[which(colnames(counts_nonref) %in% useppl),]
    counts_nonref=counts_nonref[,which(colnames(counts_nonref) %in% useppl)]
    # sum within individuals across non ref cells
    ids=unique(metadata_nonref$new_SampleID)
    new_counts_nonref=matrix(nr=nrow(counts_nonref),nc=length(ids))
    new_metadata_nonref=data.frame()
    for(o in 1:length(ids)) {
      ixi=which(metadata_nonref$new_SampleID==ids[o])
      new_counts_nonref[,o]=rowMeans(counts_nonref[,ixi,drop=FALSE]) # mean expression for each gene across all cell types
      new_metadata_nonref=rbind(new_metadata_nonref,
                                metadata_nonref %>% 
                                  filter(new_SampleID==ids[o]) %>%
                                  select(-celltype,-ncells) %>%
                                  distinct() %>%
                                  head(.,1))
    }
    counts_nonref=new_counts_nonref
    metadata_nonref=new_metadata_nonref
    # perform DEG between each cell type and all others
    # rows of metadata_ref do not match rows of metadata_nonref
    if(is.null(colnames(counts_nonref))) colnames(counts_nonref)=metadata_nonref$new_SampleID
    rownames(counts_nonref)=rownames(counts_ref)
    merged_metadata=inner_join(
      metadata_ref %>% select(new_SampleID,Age,Sex,BraakStage),
      metadata_ref %>% select(new_SampleID),by='new_SampleID')
    counts_ref=counts_ref[,merged_metadata$new_SampleID]
    counts_nonref=counts_nonref[,merged_metadata$new_SampleID]
    # loopover genes
    genes=rownames(counts_ref)
    tcounts_ref=log2(1+counts_ref)
    tcounts_nonref=log2(1+counts_nonref)
    ref_means=rowMeans(tcounts_ref)
    nonref_means=rowMeans(tcounts_nonref)
    ref_vars=apply(tcounts_ref,1,var)/ncol(counts_ref)
    nonref_vars=apply(tcounts_nonref,1,var)/ncol(counts_nonref)
    stats=(ref_means-nonref_means)^2/(ref_vars+nonref_vars)
    pvalues=pchisq(stats,1,lower.tail=FALSE)
    toadd=data.frame(
      New_Group=groups[i],
      gene=rownames(counts_ref),
      target_celltype=celltypes[j],
      mean_target=ref_means,
      mean_nontarget=nonref_means,
      var_target=ref_vars,
      var_nontarget=nonref_vars,
      statistic=stats,
      pval=pvalues)
    rdf=rbind(rdf,toadd)
  }
}
# save
setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype')
saveRDS(rdf,'DEGs_ttest.Rds')
# saveRDS(rdf,'DEGs_ttest_no_Homeostasis.Rds')
### save for Jielin
setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype')
rdf=readRDS('DEGs_ttest.Rds')
cts=unique(rdf$target_celltype)
for(i in 1:length(cts)) {
  groups=unique(rdf %>% filter(target_celltype==cts[i]) %>% pull(New_Group))
  for(j in 1:length(groups)) {
    dataout=rdf %>%
      filter(target_celltype==cts[i],New_Group==groups[j]) %>%
      mutate(logfc=log2(mean_target/mean_nontarget)) %>%
      mutate(fdr=p.adjust(pval,'BH')) %>%
      # filter(fdr<0.05) %>%
      # filter(logfc>-Inf & logfc<Inf) %>%
      # select(gene,logfc)
      select(symbol,p_val=p,log2FC=log2fc,p_val_adj=fdr)
    # setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype/for_network')
    setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/DEG_results_for_Yadi/output')
    fpout=paste0(groups[j],'_',cts[i],'_vs_others.csv')
    colnames(dataout)[1]=''
    write.table(dataout,fpout,row.names=FALSE,quote=FALSE,col.names=FALSE,sep=',')
  }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# LOCALLY!
library(dplyr);library(ggplot2)
setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/results/new/DEGs_within_group_between_microglia_subtype')
data=readRDS('DEGs_ttest.Rds')
# data=readRDS('DEGs_ttest_no_Homeostasis.Rds')
data %>%
  mutate(pval=ifelse(pval==0,min(pval[pval>0]),pval)) %>%
  mutate(mean_diff=mean_target-mean_nontarget) %>%
  mutate(fdrsig=p.adjust(pval,'BH')<0.001 & abs(mean_diff)>0.5) %>%
  ggplot(aes(x=mean_diff,y=-log10(pval),color=fdrsig)) +
  geom_point() +
  scale_color_manual(values=c('#A3D3FF','#FA9B9B')) +
  facet_grid(New_Group~target_celltype) +
  theme_bw() +
  theme(strip.background=element_rect(fill="gray95")) +
  theme(legend.text=element_text(size=11),
        legend.position='none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  labs(x='mean expression (target cell type) - mean expression (non-target cell types)',
       y=expression('-log'[10]*'(P-value)'))
setwd('/Volumes/chengflab/Cheng-Noah/DEG_testing/TACA2/plots/')
ggsave('microglia_subtypes_WithinGroupBetweenCell_DEGs.pdf',width=8,height=4)
ggsave('microglia_subtypes_WithinGroupBetweenCell_DEGs.png',width=8,height=4,dpi=400)
























