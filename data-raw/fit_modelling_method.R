#load libraries
library(QDNAseqmod)
library(Biobase)
library(dplyr)

qc.data <- read.table("data-raw/fit_QC_predownsample_withsmooth.tsv",header = T,sep = "\t")
## This data is available in the zenodo archive - https://doi.org/10.5281/zenodo.8304868
# File not provided
# Download and decompress the britroc_1_additional_data.zip and provide path to the britroc_smoothed_copyNumbersSegmented.rds file
rdsdata <- ""

if(rdsdata == ""){
  stop("no rds data provided - Download britroc_1_additional_data.zip from
  https://doi.org/10.5281/zenodo.8304868 and provide path to the britroc_smoothed_copyNumbersSegmented.rds file")
}

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

    to_use <- fData(relcn)$use
    relcn <- relcn[to_use,]

    res <- do.call(rbind,lapply(1:length(ploidies),FUN = function(i){
        ploidy <- ploidies[i]
        purity <- purities[i]

        # features calculated as in the pipeline gridsearch
        # 0.4 = default homozygous_threshold, 30k bin size (bp)
        gridVals <- rswgsabsolutecn::gridStats(obj = relcn,ploidy = ploidy,purity = purity)
        errors <- rswgsabsolutecn::getErrors(data = gridVals)
        stats <- data.frame(sample = sample,
                            segments = rswgsabsolutecn::getSegCount(relcn),
                            clonality = errors[["clonality"]],
                            homozygousLoss = sum(gridVals$abs_seg <= 0.4) * 30000,
                            segvariance = errors[["MedianSegVar"]],
                            ploidy = round(ploidy,digits = 1),
                            purity = round(purity,digits = 2),
                            use = uses[i])

        return(stats)
    }))
}))

fitstats_ord <- fitstats[order(fitstats$sample,fitstats$clonality),]

cellQC <- fitstats_ord %>% ## generated in
    tidyr::drop_na() %>%
    dplyr::mutate(use = factor(use))

## Cell QC file used for modelling
write.table(cellQC,file = "data-raw/britroc_multifit_meta.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
