library(data.table);library(dplyr);library(magrittr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions ####
pb=function(md,cm,covariates=c('Age','Sex')) {
    if(nrow(md)!=ncol(cm)) stop('number of cells in count matrix does not match number of rows in metadata')
    ids=unique(md$new_SampleID)
    pbmat=matrix(nr=nrow(cm),nc=length(ids))
    colnames(pbmat)=ids
    rownames(pbmat)=rownames(cm)
    ncells=c()
    for(i in 1:length(ids)) {
        IXI=which(md$new_SampleID==ids[i])
        cmi=as.matrix(cm[,IXI])
        pbcount=rowSums(cmi)
        pbmat[,i]=pbcount
        ncells[i]=length(IXI)
    }
    pb_metadata=as.data.frame(md)[,c('new_SampleID',covariates)] %>% 
        distinct(new_SampleID,.keep_all=TRUE) %>%
        mutate(ncells=ncells)
    # order of IDs in pb_metadata$new_SampleID will be the same as the order of IDs in `ids`
    # see: all.equal(pb_metadata$new_SampleID,ids)==TRUE
    list(pb_metadata=pb_metadata,pb_counts=pbmat)
}
consolehist=function(x,bins=10,size=10,range=c(min(x),max(x))) {
    x=x[x>=range[1] & x<=range[2]]
    cs=seq(min(x),max(x),length.out=bins+1)
    n=c()
    for(i in 2:length(cs)) {
        xi=x[x>cs[i-1] & x<=cs[i]]
        n[i-1]=length(xi)
    }
    actualn=n # when rescaling, some small counts will be rounded to 0, but I don't want to forget that they originally were not 0
    n=(n-min(n))/(max(n)-min(n))*size
    n=ceiling(n)
    prefix=c()
    for(i in 1:bins) prefix[i]=paste0(round(cs[i],2), ' (', actualn[i], ')')
    nleftchars=nchar(prefix)
    startat=max(nleftchars)+1
    res=c()
    for(i in 1:bins) {
        space=paste(rep('@',n[i]),collapse='')
        nonspace=paste(rep('-',size-n[i]),collapse='')
        toadd=paste(c(space,nonspace),collapse='')
        res[i]=toadd
    }
    # all toadd's will be of length `size` and will start at the `startat` position. need to fill from the end of `prefix`
    ntofill=startat-nleftchars
    for(i in 1:bins) {
        addedspace=paste(rep(' ',startat-nchar(prefix[i])),collapse='')
        prefix[i]=paste(prefix[i],addedspace,collapse='')
        res[i]=paste(prefix[i],res[i],collapse='')
    }
    # print
    for(i in 1:bins) cat(res[i],'\n')
}
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# the goal is to find DEGs between multiple diagnostic groups in Sst and Lamp5 cell types
setwd('/home/lorincn/isilon/Cheng-Qiu/TACA2/Integration_test/output_TACA2/Data_shared_Noah10/Sst_Neuron')
countmat=Matrix::readMM('matrix.mtx.gz')
metadata=readRDS('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/metadata/Sst_Neuron.Rds')
genes=readRDS('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/genes.Rds')
if(is.null(rownames(countmat))) rownames(countmat)=genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# remove some subjects up front
subjectix=which(metadata$new_SampleID %in% c('R6415047'))
if(length(subjectix)>0) {
    countmat=countmat[,-subjectix]
    metadata=metadata[-subjectix,]
}
# remove some datasets up front
datasetix=which(metadata$Dataset %in% c('37824655','37824663'))
if(length(datasetix)>0) {
    countmat=countmat[,-datasetix]
    metadata=metadata[-datasetix,]
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
## QC on cells and genes
cellQC=function(cm,metadata,zero_thres=0.9) {
    BOO=+(cm==0)
    BOO=colMeans(as.matrix(BOO)) # proportions of 0 read counts across genes for each cell
    BOO=which(BOO<zero_thres) # indices of cells with less than (zero_thres*100)% 0's
    out1=cm[,BOO] # keep cells with <(zero_thres*100)% zeros
    out2=metadata[BOO,] # keep cells with <(zero_thres*100)% zeros
    list(cm=out1,metadata=out2)
}
geneQC=function(countmat,zero_thres=0.9) {
    BOO=+(countmat==0)
    BOO=rowMeans(as.matrix(BOO)) # proportions of 0 read counts across cells for each gene
    BOO=which(BOO<zero_thres) # indices of genes with less than (zero_thres*100)% 0's
    countmat[BOO,]
}
# use all data - ie do not subset right now
# QC
QC=cellQC(countmat,metadata,0.9)
dat=geneQC(QC$cm,0.9)
countmat=dat
metadata=QC$metadata
# pseudobulking
covariates=c('Age','Sex','Dataset','New_Group','BraakStage')
pbdata=pb(metadata,countmat,covariates)

# clean data
countmat=pbdata$pb_counts # pseudobulk count matrix
covariatedf=pbdata$pb_metadata # pseudobulk metadata
covariatedf$Age=as.numeric(covariatedf$Age)
comparisons=list(
    ADvsCtrl=c('New_Group','AD','Normal'), # <variable containing group indicators>, <group 1 indicator>, <group 2 indicator>
    ADvsResilience=c('New_Group','AD','Resilience'),
    ADvsPART=c('New_Group','AD','PART'),
    PARTvsCtrl=c('New_Group','PART','Normal'),
    LateVsEarly=c('DiseaseStage','Early','Late')
)
for(cc in 1:length(comparisons)) {
    resultdf=data.frame()
    compvar=comparisons[[cc]][1]
    groups=comparisons[[cc]][-1]
    covariatedf$comparator=covariatedf[,compvar]
    pplix=which(covariatedf$comparator %in% groups)
    metadata_cc=covariatedf[pplix,]
    countmat_cc=countmat[,pplix]
    # x and y values
    group1ix=which(metadata_cc$comparator==groups[1])
    group2ix=which(metadata_cc$comparator==groups[2])
    x_counts=countmat_cc[,group1ix]
    x_covariates=metadata_cc[group1ix,]
    y_counts=countmat_cc[,group2ix]
    y_covariates=metadata_cc[group2ix,]
    # iterate over genes
    for(gene in 1:nrow(x_counts)) {
        thisgene=rownames(countmat)[gene]
        if(gene%%floor(nrow(x_counts)*0.01)==0) cat(round(gene/nrow(x_counts)*100),'%\n',sep='')
        # make dataframe to analyze
        dfgene=rbind(
            x_covariates %>% select(new_SampleID,Age,Sex,Dataset,ncells) %>% mutate(group=0),
            y_covariates %>% select(new_SampleID,Age,Sex,Dataset,ncells) %>% mutate(group=1)) %>%
            mutate(count=c(x_counts[gene,],y_counts[gene,])) %>%
            mutate(normcounts=count/ncells) %>%
            mutate(Age=as.numeric(Age)) %>%
            na.omit()
        # make sure analysis is possible first - ie I worry about small sample size problems
        if(var(dfgene$group,na.rm=TRUE)==0) next
        if(var(dfgene$normcounts,na.rm=T)==0) next
        # perform analysis (use mean counts for all dependent variables)
        wilcox=wilcox.test(normcounts~group,data=dfgene)
        wilcox_trans=wilcox.test(I(log2(1+normcounts))~group,data=dfgene)
        ttest=t.test(normcounts~group,data=dfgene)
        ttest_trans=t.test(I(log2(1+normcounts))~group,data=dfgene)
        linreg=lm(normcounts~Age+Sex+Dataset+group,data=dfgene)
        linreg_trans=lm(I(log2(1+normcounts))~Age+Sex+Dataset+group,data=dfgene)
        # rlinreg=MASS::rlm(normcounts~Age+Sex+Dataset+group,data=dfgene)
        # rlinreg_trans=MASS::rlm(I(log2(1+normcounts))~Age+Sex+Dataset+group,data=dfgene)
        linreg_sexdiff=lm(normcounts~Age+Sex+Dataset+group+group*Sex,data=dfgene)
        linreg_trans_sexdiff=lm(I(log2(1+normcounts))~Age+Sex+Dataset+group+group*Sex,data=dfgene)
        # R-squared
        rsq=summary(linreg)$r.squared
        rsq_trans=summary(linreg_trans)$r.squared
        # store results
        toadd=data.frame(
            gene=thisgene,
            # summary stats
            mean_group0=mean(dfgene$normcounts[dfgene$group==0],na.rm=TRUE),
            mean_group1=mean(dfgene$normcounts[dfgene$group!=0],na.rm=TRUE),
            s2bar_group0=var(dfgene$normcounts[dfgene$group==0],na.rm=TRUE),
            s2bar_group1=var(dfgene$normcounts[dfgene$group!=0],na.rm=TRUE),
            median_group0=median(dfgene$normcounts[dfgene$group==0],na.rm=TRUE),
            median_group1=median(dfgene$normcounts[dfgene$group!=0],na.rm=TRUE),
            mad_group0=mad(dfgene$normcounts[dfgene$group==0],na.rm=TRUE),
            mad_group1=mad(dfgene$normcounts[dfgene$group!=0],na.rm=TRUE),
            nppl_group0=sum(dfgene$group==0),
            nppl_group1=sum(dfgene$group!=0),
            ntotal_cells_resilience=ncol(x_counts),
            ntotal_cells_ad=ncol(y_counts),
            # Wilcoxon
            wilcox_p=wilcox$p.value,
            wilcox_trans_p=wilcox_trans$p.value,
            # t-test
            ttest_p=ttest$p.value,
            ttest_trans_p=ttest_trans$p.value,
            # linear regression
            linreg_est=tail(coef(linreg),1),
            linreg_se=tail(sqrt(diag(vcov(linreg))),1),
            linreg_trans_est=tail(coef(linreg_trans),1),
            linreg_trans_se=tail(sqrt(diag(vcov(linreg_trans))),1),
            linreg_rsquared=rsq,
            linreg_trans_rsquared=rsq_trans,
            # robust linear regression
            # rlinreg_est=tail(coef(rlinreg),1),
            # rlinreg_se=tail(sqrt(diag(vcov(rlinreg))),1),
            # rlinreg_trans_est=tail(coef(rlinreg_trans),1),
            # rlinreg_trans_se=tail(sqrt(diag(vcov(rlinreg_trans))),1)
            # tests for sex differences
            linreg_sexdiff_est=tail(coef(linreg_sexdiff),1),
            linreg_sexdiff_se=tail(sqrt(diag(vcov(linreg_sexdiff))),1),
            linreg_sexdiff_trans_est=tail(coef(linreg_trans_sexdiff),1),
            linreg_sexdiff_trans_se=tail(sqrt(diag(vcov(linreg_trans_sexdiff))),1)
        )
        resultdf=rbind(resultdf,toadd)
    }
    # save/write out results
    cat(' ',names(comparisons[cc]),'\n',sep='')
    setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/results/Sst_Neuron')
    saveRDS(resultdf,paste0(groups[1],'_vs_',groups[2],'_DEGs.Rds'))
}
