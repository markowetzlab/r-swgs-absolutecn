#load libraries
library(QDNAseqmod)
library(Biobase)
library(ggplot2)
library(stringr)
library(CINfits)

qc.data <- read.table("data-raw/fit_QC_predownsample_withsmooth.tsv",header = T,sep = "\t")
## This data is available in the zenodo archive - https://doi.org/10.5281/zenodo.8304868
## Download and decompress the britroc_1_additional_data.zip and provide path to the britroc_smoothed_copyNumbersSegmented.rds file
rdsdata <- paste0("../../../OneDrive - CRUK Cambridge Institute/britroc-cn-analysis/absolute_PRE_down_sampling/britroc_smoothed_copyNumbersSegmented.rds")

rds.rel <- readRDS(rdsdata)
samples <- qc.data[which(qc.data$SAMPLE_ID %in% colnames(rds.rel)),]

abs_profiles <- rds.rel[Biobase::fData(rds.rel)$use,]

rds.obj <- abs_profiles[,colnames(abs_profiles) %in% samples$SAMPLE_ID]

fitstats <- do.call(rbind,lapply(colnames(rds.obj),FUN = function(sample){

    ploidies <- samples$ploidy[samples$SAMPLE_ID == sample]
    purities <- samples$purity[samples$SAMPLE_ID == sample]
    # keep each fit's label with its ploidy/purity
    uses <- samples$use[samples$SAMPLE_ID == sample]

    ind<-which(colnames(rds.obj)==sample)
    relcn<-rds.obj[,ind]
    # added by PS
    to_use <- fData(relcn)$use
    relcn <- relcn[to_use,]
    # added by PS
    copynumber<-assayDataElement(relcn,"copynumber")
    rel_ploidy<-mean(copynumber,na.rm=T)

    cn<-assayDataElement(relcn,"copynumber")
    seg<-assayDataElement(relcn,"segmented")
    cn<-as.numeric(cn[!is.na(cn),])
    seg<-as.numeric(seg[!is.na(seg),])

    res <- do.call(rbind,lapply(1:length(ploidies),FUN = function(i){
        ploidy <- ploidies[i]
        purity <- purities[i]

        cellploidy<-ploidy*purity+2*(1-purity)
        seqdepth<-rel_ploidy/cellploidy

        abs_seg <- depthtocn(seg,purity,seqdepth)
        abs_cnbin <- depthtocn(cn,purity,seqdepth)

        segmentVar <- calculateSegmentVar(abs_seg,abs_cnbin)

        relcn_temp <- relcn
        assayDataElement(relcn_temp,"segmented")[,1] <- abs_seg

        tab <- getSegTable(relcn_temp)

        stats <- as.data.frame(cbind(CINfits::calculateCINStats(tab),segvariance=segmentVar))

        stats <- stats %>% tibble::rownames_to_column("sample") %>%
            dplyr::mutate(grid_ploidy = round(rep(ploidy,times=nrow(.)),digits = 1)) %>%
            dplyr::mutate(grid_purity = round(rep(purity,times=nrow(.)),digits = 2)) %>%
            dplyr::mutate(use = uses[i])

        return(stats)
    }))
}))

fitstats_ord <- fitstats[order(fitstats$sample,fitstats$clonality),]

cellQC <- fitstats_ord %>% ## generated in
    tidyr::drop_na() %>%
    dplyr::select(-ploidy) %>%
    dplyr::rename(ploidy = grid_ploidy) %>%
    dplyr::rename(purity = grid_purity) %>%
    dplyr::mutate(use = factor(use))

## Cell QC file used for modelling
write.table(cellQC,file = "data-raw/britroc_multifit_meta.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
