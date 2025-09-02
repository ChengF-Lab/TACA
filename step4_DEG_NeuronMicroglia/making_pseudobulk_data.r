library(data.table);library(dplyr);library(magrittr)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# functions ####
pb=function(md,cm,covariates=c('Age','Sex'),subtype_column=NULL) {
    if(nrow(md)!=ncol(cm)) stop('number of cells in count matrix does not match number of rows in metadata')
    if(is.null(subtype_column)) {
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
    } else {
        pbmat=matrix(nr=nrow(cm),nc=1) # will overwrite this column later
        rownames(pbmat)=rownames(cm)
        rmeta=data.frame()
        sts=unique(md %>% pull(subtype_column))
        k=0;nn=c()
        for(j in 1:length(sts)) {
            ix=metadata[,subtype_column]==sts[j]
            ix=which(c(ix))
            mdj=metadata[ix,]
            cmj=cm[,ix,drop=FALSE]
            ids=unique(mdj$new_SampleID)
            for(i in 1:length(ids)) {
                k=k+1
                nn[k]=paste0(sts[j],'_',ids[i])
                IXI=which(mdj$new_SampleID==ids[i])
                cmi=as.matrix(cmj[,IXI])
                pbcount=rowSums(cmi)
                if(k==1) pbmat[,k]=pbcount else pbmat=cbind(pbmat,pbcount)
                toadd=as.data.frame(mdj)[,c('new_SampleID',covariates)] %>% 
                    filter(new_SampleID==ids[i]) %>%
                    distinct(new_SampleID,.keep_all=TRUE) %>%
                    mutate(ncells=length(IXI),subtype=sts[j])
                rmeta=rbind(rmeta,toadd)
            }
            colnames(pbmat)=nn
            pb_metadata=rmeta
        }
    }
    # order of IDs in pb_metadata$new_SampleID will be the same as the order of IDs in `ids`
    # see: all.equal(pb_metadata$new_SampleID,ids)==TRUE
    list(pb_metadata=pb_metadata,pb_counts=pbmat)
}
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
# dictionary of microglia subtypes
microglia_type_list=list(
    DAM=c('DAM1','DAM2'),
    IFN='IFN Microglia',
    Inflammation='Inflammation Microglia',
    Tau='Tau Microglia',
    TCell='Leukocyte_T cell',
    Monocyte='Monocyte',
    Macrophage='Macrophage',
    Proliferation='Proliferation Microglia',
    MHC=c('MHC Microglia1','MHC Microglia2'),
    Homeostasis=c('Homeostasis Microglia1','Homeostasis Microglia2')
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# simply save pseudobulk data the user wants ####
celltype='Neuron0' # Astrocyte,Immune,Oligodedrocyte,OPC,Others,Sst_Neuron,Neuron0,Neuron1,...,Neuron23
microgliasubtype='none' # DAM,IFN,Inflammation,Tau,TCell,Macrophage,Proliferation, MHC
usecohorts='all' # all, none
dropcohorts='37824655,37824663' # all, none
removesubjects='R6415047' # <new_SampleID>,NA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
setwd('/home/lorincn/isilon/Cheng-Qiu/TACA2/Integration_test/output_TACA2/Data_shared_Noah11') # setwd('/home/lorincn/isilon/Cheng-Qiu/TACA2/Integration_test/output_TACA2/Data_shared_Noah10')
if(!(celltype %in% dir())) stop(paste0('celltype ',celltype,' not found in ',getwd()))
setwd(celltype)
countmat=Matrix::readMM('matrix.mtx.gz')
metafp=paste0('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/metadata/',celltype,'.Rds')
metadata=readRDS(metafp)
genes=readRDS('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/genes.Rds')
if(is.null(rownames(countmat))) rownames(countmat)=genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# filter to only certain datasets?
usecohorts=unlist(strsplit(usecohorts,','))
dropcohorts=unlist(strsplit(dropcohorts,','))
if(!all(usecohorts=='all')) {
    datasetix=which(metadata$Dataset %in% usecohorts)
    countmat=countmat[,datasetix]
    metadata=metadata[datasetix,]
}
if(!all(dropcohorts=='none')) {
    datasetix=which(metadata$Dataset %in% dropcohorts)
    if(length(datasetix)>0) {
        countmat=countmat[,-datasetix]
        metadata=metadata[-datasetix,]
    }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# remove some subjects?
if(!all(removesubjects=='none')) {
    removesubjects=unlist(strsplit(removesubjects,','))
    subjectix=which(metadata$new_SampleID %in% removesubjects)
    if(length(subjectix)>0) {
        countmat=countmat[,-subjectix]
        metadata=metadata[-subjectix,]
    }
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# filter to just certain cell subtypes?
if(!all(microgliasubtype=='none') & celltype!='Immune') stop('cannot specify microglia cell subtypes for non-immune cells')
if(!all(microgliasubtype=='none')) {
    # pseudobulk for each cell type subtype
    microgliasubtype=unlist(strsplit(microgliasubtype,','))
    # make sure all are in data
    # indata=c(); for(i in 1:length(microgliasubtype)) indata[i]='MHC Microglia1' %in% unique(metadata$L3_label) # just specific L3 subtypes (for running manually)
    indata=c(); for(i in 1:length(microgliasubtype)) indata[i]=any(microglia_type_list[[i]] %in% unique(metadata$L3_label))
    microgliasubtype=microgliasubtype[which(indata)]
    typeix=which(names(microglia_type_list) %in% microgliasubtype)
    actual_types=list(); for(j in 1:length(typeix)) actual_types[[j]]=microglia_type_list[[typeix[j]]]
    names(actual_types)=microgliasubtype
    for(i in 1:length(actual_types)) {
        ixi=which(metadata$L3_label %in% actual_types[[i]])
        countmati=countmat[,ixi]
        metadatai=metadata[ixi,]
        ## cell and gene QC
        QC=cellQC(countmati,metadatai,0.99)
        metadatai=QC$metadata
        QC=geneQC(QC$cm,0.99)
        countmati=QC
        # actual pseudobulking
        covariates=c('Age','Sex','Dataset','New_Group','BraakStage','DiseaseStage')
        pbdata=pb(metadatai,countmati,covariates)
        # save
        setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data')
        setwd(celltype)
        fpout=paste(usecohorts,collapse='_')
        fpout=paste0(fpout,'_cohorts')
        fpout=paste0(fpout,'_',microgliasubtype[i],'_subtype')
        if(!(all(dropcohorts=='none'))) fpout=paste0(fpout,'_',paste(dropcohorts,collapse='_'),'_droppedCohorts')
        if(!(all(removesubjects=='none'))) fpout=paste0(fpout,'_',paste(removesubjects,collapse='_'),'_droppedSubjects')
        fpout=paste0(fpout,'_99pctQC.Rds')
        saveRDS(pbdata,fpout)
    }
} else {
    # pseudobulk across all cell type subtypes
    ## cell and gene QC
    QC=cellQC(countmat,metadata,0.99)
    metadata=QC$metadata
    QC=geneQC(QC$cm,0.99)
    countmat=QC
    # actual pseudobulking
    covariates=c('Age','Sex','Dataset','New_Group','BraakStage','DiseaseStage')
    # pbdata=pb(metadata,countmat,covariates,subtype_column='L2_label') # L2_label for Neuron0-23
    pbdata=pb(metadata,countmat,covariates) # L2_label for Neuron0-23
    # save
    setwd('/home/lorincn/isilon/Cheng-Noah/DEG_testing/TACA2/pseudobulk_data')
    setwd(celltype)
    fpout=paste(usecohorts,collapse='_')
    fpout=paste0(fpout,'_cohorts')
    if(!(all(dropcohorts=='none'))) fpout=paste0(fpout,'_',paste(dropcohorts,collapse='_'),'_droppedCohorts')
    if(!(all(removesubjects=='none'))) fpout=paste0(fpout,'_',paste(removesubjects,collapse='_'),'_droppedSubjects')
    fpout=paste0(fpout,'_99pct_QC.Rds')
    saveRDS(pbdata,fpout)
}

