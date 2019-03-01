###### Phosphosites #####

setwd("H:/Documents/05_collaborations_small_projects/Umarah_phosphoproteomics_normalization")
setwd("/home/thomas/Documents/PhD/collaborations/phosphoproteomics_umarah/")
x <- read.table("valid_filt_PhosphoProteomics_0_64_copy.txt", header=TRUE, fill=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(x) <- paste(x[,"Protein"], x[,"Position"], x[,"Amino.acid"], rownames(x), sep="_")
x_prot <- read.table("valid_filt_Proteomics_0_64_copy.txt", header=TRUE, fill=TRUE, sep="\t", stringsAsFactors=FALSE)
x_prot <- x_prot[-1, ] ## remove first row

## take first entry of "Majority.protein.IDs"
rownames(x_prot) <- unlist(lapply(strsplit(x_prot[, "Majority.protein.IDs"], split=";"), "[", 1))

## rename ..._0_5:... to ..._05_...
colnames(x) <- gsub(x=colnames(x), pattern="[.]C_0_5", replacement=".C_05") 
colnames(x) <- gsub(x=colnames(x), pattern="[.]R_0_5", replacement=".R_05")
colnames(x_prot) <- gsub(x=colnames(x_prot), pattern="[.]C_0_5", replacement=".C_05") 
colnames(x_prot) <- gsub(x=colnames(x_prot), pattern="[.]R_0_5", replacement=".R_05") 

## rename 
colnames(x)[colnames(x) == "Intensity.R_16_22"] <- "Intensity.R_16_2_2"

x <- x[,grep(colnames(x), pattern="Intensity")]
x_prot <- x_prot[, grep(colnames(x_prot), pattern="Intensity")]

## remove specific colnames 
x <- x[, !colnames(x) %in% c("Intensity",  "Intensity.0_6", "Intensity.0_7")]
x_prot <- x_prot[, -c(which(colnames(x_prot)=="Intensity"),  grep(colnames(x_prot), pattern="y[.]0_6|y[.]0_7|y[.]0_8"))]

x <- as.matrix(x)
mode(x) <- "numeric"
x_prot <- as.matrix(x_prot)
mode(x_prot) <- "numeric"

# plot(hclust(dist(x_prot_imp)), labels=F)
# x_prot_hclust <- hclust(dist(x_prot_imp))
# cutree_3 <- cutree(x_prot_hclust, k=12)
# which(cutree_3 != 1)
## manually remove outliers in proteome file
outlier <- read.table("outlier features_proteome.csv", sep=";", header=TRUE, stringsAsFactors=FALSE)[, "Name"]
na_outlier <- apply(x_prot[which(rownames(x_prot) %in% outlier),], 1, function(x) sum(is.na(x)))
hist(na_outlier)
x_prot <- x_prot[-which(rownames(x_prot) %in% outlier),]

## divide by 75% quantile and 50% quantile



ind_batch1 <- grep(colnames(x_prot), pattern="0_1|C_05_1|C_1_1|C_2_1|C_4_1|C_8_1|C_16_1|C_32_1|C_64_1|R_05_1|R_1_1|R_2_1|R_4_1|R_8_1|R_16_1|R_32_1|R_64_1")
ind_batch2 <- grep(colnames(x_prot), pattern="0_2|C_05_2|C_1_2|C_2_2|C_4_2|C_8_2|C_16_2|C_32_2|C_64_2|R_05_2|R_1_2|R_2_2|R_4_2|R_8_2|R_16_2|R_32_2|R_64_2")
ind_batch3 <- grep(colnames(x_prot), pattern="0_3|C_05_3|C_1_3|C_2_3|C_4_3|C_8_3|C_16_3|C_32_3|C_64_3|R_05_3|R_1_3|R_2_3|R_4_3|R_8_3|R_16_3|R_32_3|R_64_3")
ind_batch4 <- grep(colnames(x_prot), pattern="0_4|C_05_4|C_1_4|C_2_4|C_4_4|C_8_4|C_16_4|C_32_4|C_64_4|R_05_4|R_1_4|R_2_4|R_4_4|R_8_4|R_16_4|R_32_4|R_64_4")
ind_batch5 <- grep(colnames(x_prot), pattern="0_5|C_05_5|C_1_5|C_2_5|C_4_5|C_8_5|C_16_5|C_32_5|C_64_5|R_05_5|R_1_5|R_2_5|R_4_5|R_8_5|R_16_5|R_32_5|R_64_5")
ind_batch_l <- c(ind_batch1, ind_batch2, ind_batch3, ind_batch4, ind_batch5)
batch_prot <- character(length(unlist(ind_batch_l)))
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.0_1"))] <- "red"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.?_*_1*"))] <- "red"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.0_2*"))] <- "blue"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.?_*_2*"))] <- "blue"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.0_3*"))] <- "green"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.?_*_3*"))] <- "green"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.0_4*"))] <- "yellow"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.?_*_4*"))] <- "yellow"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.0_5*"))] <- "orange"
batch_prot[grep(colnames(x_prot), pattern=glob2rx("Intensity.?_*_5*"))] <- "orange"




## only use replicates 1-3 for phospoproteomics
x <- x[, -grep(colnames(x), pattern=glob2rx("Intensity.?_*_4*"))]
x <- x[, -grep(colnames(x), pattern=glob2rx("Intensity.?_*_5*"))]
x <- x[, -grep(colnames(x), pattern=glob2rx("Intensity.0_4*"))]
x <- x[, -grep(colnames(x), pattern=glob2rx("Intensity.0_5*"))]

##x <- apply(x, 2, function(y) y / quantile(y, .75, na.rm=TRUE))
##x_prot <- apply(x_prot, 2, function(y) y / quantile(y, .5, na.rm=TRUE))
ind_batch1 <- grep(colnames(x), pattern="0_1|C_05_1|C_1_1|C_2_1|C_4_1|C_8_1|C_16_1|C_32_1|C_64_1|R_05_1|R_1_1|R_2_1|R_4_1|R_8_1|R_16_1|R_32_1|R_64_1")
ind_batch2 <- grep(colnames(x), pattern="0_2|C_05_2|C_1_2|C_2_2|C_4_2|C_8_2|C_16_2|C_32_2|C_64_2|R_05_2|R_1_2|R_2_2|R_4_2|R_8_2|R_16_2|R_32_2|R_64_2")
ind_batch3 <- grep(colnames(x), pattern="0_3|C_05_3|C_1_3|C_2_3|C_4_3|C_8_3|C_16_3|C_32_3|C_64_3|R_05_3|R_1_3|R_2_3|R_4_3|R_8_3|R_16_3|R_32_3|R_64_3")
##ind_batch4 <- grep(colnames(x), pattern="0_4|05_4|1_4|2_4|4_4|8_4|16_4|32_4|64_4")
##ind_batch5 <- grep(colnames(x), pattern="0_5|05_5|1_5|2_5|4_5|8_5|16_5|32_5|64_5")
ind_batch_l <- c(ind_batch1, ind_batch2, ind_batch3)##, ind_batch4, ind_batch5)
batch <- character(length(unlist(ind_batch_l)))
batch[grep(colnames(x), pattern=glob2rx("Intensity.0_1*"))] <- "red"
batch[grep(colnames(x), pattern=glob2rx("Intensity.?_*_1*"))] <- "red"
batch[grep(colnames(x), pattern=glob2rx("Intensity.0_2*"))] <- "blue"
batch[grep(colnames(x), pattern=glob2rx("Intensity.?_*_2*"))] <- "blue"
batch[grep(colnames(x), pattern=glob2rx("Intensity.0_3*"))] <- "green"
batch[grep(colnames(x), pattern=glob2rx("Intensity.?_*_3*"))] <- "green"


normalization <- "quantile" ## quantile

if (normalization == "batchEffect") {
## removeBatchEffect
    x_prot <- limma::removeBatchEffect(x_prot, batch=batch_prot)
    x <- limma::removeBatchEffect(x, batch=batch)

}

if (normalization == "quantile") {
    ## 50% and 75% quantile normalization
    x <- apply(x, 2, function(y) y / quantile(y, .75, na.rm=TRUE))
    x_prot <- apply(x_prot, 2, function(y) y / quantile(y, .5, na.rm=TRUE))    
}




par(mfrow=c(5,5))
for(i in 1:ncol(x)) hist(log2(x[,i]))
for(i in 1:ncol(x_prot)) hist(log2(x_prot[,i]))



## check for outliers --> PCA and manually remove
library(pcaMethods)
x_outl <- x
x_prot_outl <- x_prot
x_outl <- ifelse(is.na(x_outl), NA, x)
x_prot_outl <- ifelse(is.na(x_prot_outl), NA, x_prot)
rpcaNN <- pca(t(x_outl), method="ppca", center=TRUE, scale="none", nPcs=5, seed=3455, na.rm=TRUE)
rpcaNN_prot <- pca(t(x_prot_outl), method="ppca", center=TRUE, scale="none", nPcs=5, seed=3455, na.rm=TRUE)
par(mfrow=c(1,1))
slplot(rpcaNN, pcs=c(1,2), scoresLoadings=c(T, F), scex=0.5) ## looks ok, do not remove outliers
slplot(rpcaNN_prot, pcs=c(1,2), scoresLoadings=c(T, F), scex=0.5) 


## filtering by timepoints, if 3/5 we take it
get_NA <- function(x, pattern="Intensity[.]0_") {
    inds <- grep(colnames(x), pattern=pattern) ## 1|2, 3|4, 5|6, 7, 8   
    if (pattern == "Intensity[.]0_") {
        inds_rep <- unlist(lapply(strsplit(colnames(x)[inds], split="_"), "[", 2))
        inds_rep <- lapply(strsplit(inds_rep, split="[.]"), "[", 1)
    } else {
        inds_rep <- unlist(lapply(strsplit(colnames(x)[inds], split="_"), "[", 3))   
        inds_rep <- lapply(strsplit(inds_rep, split="[.]"), "[", 1)
    }
    
    get_NA_rep <- function(rep) {
        if (sum(inds_rep==rep) == 1) {
            return(!is.na(x[,inds[inds_rep == rep]]))
        } else {
            return(!apply(x[,inds[inds_rep == rep]], 1, function(x) all(is.na(x))))
        }    
    }
    t_1 <- get_NA_rep(rep="1")
    t_2 <- get_NA_rep(rep="2")
    t_3 <- get_NA_rep(rep="3")
    t_4 <- get_NA_rep(rep="4")
    t_5 <- get_NA_rep(rep="5")
    threshold <- 3
    if (length(unique(inds_rep)) != 5) {print("set threshold to 2"); threshold <- 2}
    t <- t_1 + t_2 + t_3 + t_4 + t_5 >= threshold
    return(t)
}

t0 <- get_NA(x=x, pattern="Intensity[.]0_")
tC05 <- get_NA(x=x, pattern="Intensity[.]C_05_") 
tR05 <- get_NA(x=x, pattern="Intensity[.]R_05_")
tC1 <- get_NA(x=x, pattern="Intensity[.]C_1_") 
tR1 <- get_NA(x=x, pattern="Intensity[.]R_1_") 
tC2 <- get_NA(x=x, pattern="Intensity[.]C_2_") 
tR2 <- get_NA(x=x, pattern="Intensity[.]R_2_") 
tC4 <- get_NA(x=x, pattern="Intensity[.]C_4_") 
tR4 <- get_NA(x=x, pattern="Intensity[.]R_4_") 
tC8 <- get_NA(x=x, pattern="Intensity[.]C_8_")
tR8 <- get_NA(x=x, pattern="Intensity[.]R_8_") ## no R_8_4 
tC16 <- get_NA(x=x, pattern="Intensity[.]C_16_")
tR16 <- get_NA(x=x, pattern="Intensity[.]R_16_")
tC32 <- get_NA(x=x, pattern="Intensity[.]C_32_") ## no C_32_5
tR32 <- get_NA(x=x, pattern="Intensity[.]R_32_")
tC64 <- get_NA(x=x, pattern="Intensity[.]C_64_") 
tR64 <- get_NA(x=x, pattern="Intensity[.]R_64_")

t0_prot <- get_NA(x=x_prot, pattern="Intensity[.]0_")
tC05_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_05_") 
tR05_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_05_")
tC1_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_1_") 
tR1_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_1_") 
tC2_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_2_") ## no C_2_5
tR2_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_2_") ## no R_2_3, R_2_4
tC4_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_4_") 
tR4_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_4_") ## no R_4_1
tC8_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_8_")
tR8_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_8_") 
tC16_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_16_")
tR16_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_16_")
tC32_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_32_") 
tR32_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_32_")
tC64_prot <- get_NA(x=x_prot, pattern="Intensity[.]C_64_") 
tR64_prot <- get_NA(x=x_prot, pattern="Intensity[.]R_64_")



## take these lines that are present in at least one of the timepoints with >= 3
t_keep <- t0 + tC05 + tR05 + tC1 + tR1 + tC2 + tR2 + tC4 + tR4 + tC8 + tR8 + tC16 + tR16 + tC32 + tR32 + tC64 + tR64 > 0
t_keep_prot <- t0_prot + tC05_prot + tR05_prot + tC1_prot + tR1_prot + tC2_prot + tR2_prot + tC4_prot + tR4_prot + tC8_prot + tR8_prot + tC16_prot + tR16_prot + tC32_prot + tR32_prot + tC64_prot + tR64_prot > 0
x <- x[t_keep, ]
x_prot <- x_prot[t_keep_prot, ]
x_outl <- x
x_outl_prot <- x_prot
x_outl <- ifelse(is.na(x_outl), NA, x)
x_outl_prot <- ifelse(is.na(x_outl_prot), NA, x_prot)
rpcaNN <- pca(t(x_outl), method="ppca", center=TRUE, scale="none", nPcs=5, seed=3455, na.rm=TRUE)
rpcaNN_prot <- pca(t(x_outl_prot), method="ppca", center=TRUE, scale="none", nPcs=5, seed=3455, na.rm=TRUE)
slplot(rpcaNN, pcs=c(1,2), scoresLoadings=c(T, F), scex=0.5) ## looks ok, do not remove outliers
slplot(rpcaNN_prot, pcs=c(1,2), scoresLoadings=c(T, F), scex=0.5) ## looks ok, do not remove outliers


## calculating average of replicates and delete unnecessary replicates
average_rep <- function(x=x, pattern_grep="Intensity[.]0_1", pattern_bind="Intensity.0_1") {
    tmp <- apply(x[, grep(colnames(x), pattern=pattern_grep)], 1, mean, na.rm=TRUE)
    x <- x[,-grep(colnames(x), pattern=pattern_grep)]; x <- cbind(x, "new"=tmp)
    colnames(x)[colnames(x) == "new"] <- pattern_bind
    return(x)
}
x <- average_rep(x=x, pattern_grep="Intensity[.]0_1", pattern_bind="Intensity.0_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]0_2", pattern_bind="Intensity.0_2")
x <- average_rep(x=x, pattern_grep="Intensity[.]0_3", pattern_bind="Intensity.0_3")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_05_1", pattern_bind="Intensity.C_05_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_05_2", pattern_bind="Intensity.C_05_2")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_1_1", pattern_bind="Intensity.C_1_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_16_1", pattern_bind="Intensity.C_16_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_2_1", pattern_bind="Intensity.C_2_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_2_2", pattern_bind="Intensity.C_2_2")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_32_1", pattern_bind="Intensity.C_32_1")
##x <- average_rep(x=x, pattern_grep="Intensity[.]C_32_4", pattern_bind="Intensity.C_32_4")
x <- average_rep(x=x, pattern_grep="Intensity[.]C_4_1", pattern_bind="Intensity.C_4_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_05_1", pattern_bind="Intensity.R_05_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_1_1", pattern_bind="Intensity.R_1_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_16_2", pattern_bind="Intensity.R_16_2")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_2_1", pattern_bind="Intensity.R_2_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_2_2", pattern_bind="Intensity.R_2_2")
##x <- average_rep(x=x, pattern_grep="Intensity[.]R_32_5", pattern_bind="Intensity.R_32_5")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_4_1", pattern_bind="Intensity.R_4_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_64_1", pattern_bind="Intensity.R_64_1")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_64_2", pattern_bind="Intensity.R_64_2")
x <- average_rep(x=x, pattern_grep="Intensity[.]R_8_1", pattern_bind="Intensity.R_8_1")

colnames(x)[which(colnames(x) == "Intensity.C_4_2_180225173818")] <-  "Intensity.C_4_2"
x <- x[, !colnames(x) %in% "Intensity.0_2_180331041935"]
##x <- x[, !colnames(x) %in% "Intensity.C_8_4"] ## remove since there is no Intensity.R_8_4
##x <- x[, !colnames(x) %in% "Intensity.R_32_5"] ## remove since there is no Intensity.C_32_5
x <- x[, sort(colnames(x))]

x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]0_3", pattern_bind="Intensity.0_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]0_4", pattern_bind="Intensity.0_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_05_2", pattern_bind="Intensity.C_05_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_05_3", pattern_bind="Intensity.C_05_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_05_4", pattern_bind="Intensity.C_05_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_05_5", pattern_bind="Intensity.C_05_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_1_2", pattern_bind="Intensity.C_1_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_1_4", pattern_bind="Intensity.C_1_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_1_5", pattern_bind="Intensity.C_1_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_16_1", pattern_bind="Intensity.C_16_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_16_2", pattern_bind="Intensity.C_16_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_16_3", pattern_bind="Intensity.C_16_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_2_1", pattern_bind="Intensity.C_2_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_2_2", pattern_bind="Intensity.C_2_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_2_3", pattern_bind="Intensity.C_2_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_2_4", pattern_bind="Intensity.C_2_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_32_1", pattern_bind="Intensity.C_32_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_32_3", pattern_bind="Intensity.C_32_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_4_2", pattern_bind="Intensity.C_4_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_4_3", pattern_bind="Intensity.C_4_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_4_4", pattern_bind="Intensity.C_4_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_64_1", pattern_bind="Intensity.C_64_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_64_3", pattern_bind="Intensity.C_64_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_64_4", pattern_bind="Intensity.C_64_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_64_5", pattern_bind="Intensity.C_64_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]C_8_1", pattern_bind="Intensity.C_8_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_05_1", pattern_bind="Intensity.R_05_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_05_2", pattern_bind="Intensity.R_05_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_05_5", pattern_bind="Intensity.R_05_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_1_2", pattern_bind="Intensity.R_1_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_1_3", pattern_bind="Intensity.R_1_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_1_5", pattern_bind="Intensity.R_1_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_16_1", pattern_bind="Intensity.R_16_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_16_3", pattern_bind="Intensity.R_16_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_2_5", pattern_bind="Intensity.R_2_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_32_2", pattern_bind="Intensity.R_32_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_32_3", pattern_bind="Intensity.R_32_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_32_4", pattern_bind="Intensity.R_32_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_32_5", pattern_bind="Intensity.R_32_5")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_4_2", pattern_bind="Intensity.R_4_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_4_3", pattern_bind="Intensity.R_4_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_64_1", pattern_bind="Intensity.R_64_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_64_3", pattern_bind="Intensity.R_64_3")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_64_4", pattern_bind="Intensity.R_64_4")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_8_1", pattern_bind="Intensity.R_8_1")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_8_2", pattern_bind="Intensity.R_8_2")
x_prot <- average_rep(x=x_prot, pattern_grep="Intensity[.]R_8_3", pattern_bind="Intensity.R_8_3")

colnames(x_prot)[which(colnames(x_prot) == "Intensity.C_05_1_repeat")] <-  "Intensity.C_05_1"
x_prot <- x_prot[, !colnames(x_prot) %in% "Intensity.0_2_180331041935"] ## remove since there is no Intensity.R_1_4
#x_prot <- x_prot[, !colnames(x_prot) %in% "Intensity.R_2_5"] ## remove since there is no Intensity.C_2_5



## save x in x_old, x_prot in x_prot_old
x_old <- x 
x_prot_old <- x_prot

## define functions that are used later
## impute values for missing values --> calculate the minimum value per row
impute_min <- function(x) {
    x_imp <- x
    for (i in 1:nrow(x_imp)) {
        if (any(is.na(x[i,]))) {
            x_imp[i, is.na(x[i,])] <- min(x[i,], na.rm=T)
        }
        ## if all elements are equal add random noise
        if (all(x_imp[i,] - mean(x_imp[i,]) == 0)) {
             x_imp[i,] <- x_imp[i,] + abs(rnorm(ncol(x_imp), mean=0, sd=1e-10))
        }
    }
    return(x_imp)
}
## calculate_fc is used in the strategies A
calculate_fc <- function(x, pattern="^Intensity.*_05_") {
    ind <- grep(colnames(x), pattern=pattern)
    beg <- seq(1, length(ind) / 2) 
    end <- beg + length(ind) / 2
    fc <- x[, ind[end]] / x[, ind[beg]]
    return(fc)
} 
## average_cond is used in the strategies B
average_cond <- function(x, condition) {
    x_cond <- x[, grep(colnames(x), pattern=condition)]
    x_cond_mean <- apply(x_cond, 1, mean, na.rm=TRUE)
    return(x_cond_mean)
}



## Strategy 1 (only 75%/50% quantile + log2)
######### 1A  ##########
## fc per replicate
x_old <- x 
x_prot_old <- x_prot
x75 <- x <- log2(x) ## +1: do this since it can be that values are 0
x_prot50 <- x_prot <- log2(x_prot)

par(mfrow=c(5,5))
for(i in 1:ncol(x75)) hist(x75[,i])
for(i in 1:ncol(x_prot50)) hist(x_prot50[,i])

## missing samples: C_2_5, R_2_3, R_2_4, R_1_4, R_4_1, 
## impute based on mean of the remaining in that group
x_prot50 <- cbind(x_prot50, "Intensity.C_2_5"=apply(x_prot50[, c("Intensity.C_2_1", "Intensity.C_2_2", "Intensity.C_2_3", "Intensity.C_2_4")], 1, mean, na.rm=TRUE))
x_prot50 <- cbind(x_prot50, "Intensity.R_2_3"=apply(x_prot50[, c("Intensity.R_2_1", "Intensity.R_2_2", "Intensity.R_2_5")], 1, mean, na.rm=TRUE))
x_prot50 <- cbind(x_prot50, "Intensity.R_2_4"=apply(x_prot50[, c("Intensity.R_2_1", "Intensity.R_2_2", "Intensity.R_2_5")], 1, mean, na.rm=TRUE))
x_prot50 <- cbind(x_prot50, "Intensity.R_1_4"=apply(x_prot50[, c("Intensity.R_1_1", "Intensity.R_1_2", "Intensity.R_1_3", "Intensity.R_1_5")], 1, mean, na.rm=TRUE))
x_prot50 <- cbind(x_prot50, "Intensity.R_4_1"=apply(x_prot50[, c("Intensity.R_4_2", "Intensity.R_4_3", "Intensity.R_4_5", "Intensity.R_4_5")], 1, mean, na.rm=TRUE))
x_prot50 <- x_prot50[, sort(colnames(x_prot50))]




## impute missing values (use min per row)
x75_imp <- impute_min(x75)
x_prot50_imp <- impute_min(x_prot50)


# shapiro_filter <- function(x) {
#     ## conduct shapiro-wilk test, H0: underlying normal distribution
#     shapiro <- lapply(1:nrow(x), function(y) shapiro.test(x[y,]))
#     shapiro <- lapply(1:length(shapiro), function(y) shapiro[[y]][["p.value"]])
#     shapiro <- unlist(shapiro)
#     shapiro <- p.adjust(shapiro, method = "BH")
#     ## return those that are non-significant (i.e. those that are normally distributed)
#     x <- x[shapiro > 0.05,]
#     return(x)
# }
# x75_imp <- shapiro_filter(x75_imp)
# x_prot50_imp <- shapiro_filter(x_prot50_imp)



##write.table(x_prot, file="prot_median_log_normalized.csv", sep="\t", dec=".")

values0 <- x75_imp[, grep(colnames(x75_imp), pattern="Intensity.0_")]
fc_0_a1 <- values0 / values0
fc_05_a1 <- calculate_fc(x75_imp, "^Intensity.*_05_")
fc_1_a1 <- calculate_fc(x75_imp, "^Intensity.*_1_") 
fc_2_a1 <- calculate_fc(x75_imp, "^Intensity.*_2_")
fc_4_a1 <- calculate_fc(x75_imp, "^Intensity.*_4_")
fc_8_a1 <- calculate_fc(x75_imp, "^Intensity.*_8_")
fc_16_a1 <- calculate_fc(x75_imp, "^Intensity.*_16_")
fc_32_a1 <- calculate_fc(x75_imp, "^Intensity.*_32_")
fc_64_a1 <- calculate_fc(x75_imp, "^Intensity.*_64_")

fc_0_a1_m <- apply(fc_0_a1, 1, mean, na.rm=TRUE)
fc_05_a1_m <- apply(fc_05_a1, 1, mean, na.rm=TRUE)
fc_1_a1_m <- apply(fc_1_a1, 1, mean, na.rm=TRUE)
fc_2_a1_m <- apply(fc_2_a1, 1, mean, na.rm=TRUE)
fc_4_a1_m <- apply(fc_4_a1, 1, mean, na.rm=TRUE)
fc_8_a1_m <- apply(fc_8_a1, 1, mean, na.rm=TRUE)
fc_16_a1_m <- apply(fc_16_a1, 1, mean, na.rm=TRUE)
fc_32_a1_m <- apply(fc_32_a1, 1, mean, na.rm=TRUE)
fc_64_a1_m <- apply(fc_64_a1, 1, mean, na.rm=TRUE)
fc_all_a1 <- cbind(fc_0_a1_m, fc_05_a1_m, fc_1_a1_m, fc_2_a1_m, fc_4_a1_m, fc_8_a1_m, fc_16_a1_m, fc_32_a1_m, fc_64_a1_m)

values0_prot <- x_prot50_imp[, grep(colnames(x_prot50_imp), pattern="Intensity.0_")]
fc_0_a1_prot <- values0_prot / values0_prot
fc_05_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_05_")
fc_1_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_1_")
fc_2_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_2_")
fc_4_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_4_")
fc_8_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_8_")
fc_16_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_16_")
fc_32_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_32_")
fc_64_a1_prot <- calculate_fc(x_prot50_imp, "^Intensity.*_64_")

fc_0_a1_m_prot <- apply(fc_0_a1_prot, 1, mean, na.rm=TRUE)
fc_05_a1_m_prot <- apply(fc_05_a1_prot, 1, mean, na.rm=TRUE)
fc_1_a1_m_prot <- apply(fc_1_a1_prot, 1, mean, na.rm=TRUE)
fc_2_a1_m_prot <- apply(fc_2_a1_prot, 1, mean, na.rm=TRUE)
fc_4_a1_m_prot <- apply(fc_4_a1_prot, 1, mean, na.rm=TRUE)
fc_8_a1_m_prot <- apply(fc_8_a1_prot, 1, mean, na.rm=TRUE)
fc_16_a1_m_prot <- apply(fc_16_a1_prot, 1, mean, na.rm=TRUE)
fc_32_a1_m_prot <- apply(fc_32_a1_prot, 1, mean, na.rm=TRUE)
fc_64_a1_m_prot <- apply(fc_64_a1_prot, 1, mean, na.rm=TRUE)
fc_all_a1_prot <- cbind(fc_0_a1_m_prot, fc_05_a1_m_prot, fc_1_a1_m_prot, fc_2_a1_m_prot, fc_4_a1_m_prot, fc_8_a1_m_prot, fc_16_a1_m_prot, fc_32_a1_m_prot, fc_64_a1_m_prot)

##### Strategy 1 B-1 average first replicates then calculate fc ######
fc_0_b1 <- average_cond(x75_imp, "Intensity[.]0_") / average_cond(x75_imp, "Intensity[.]0_")
fc_05_b1 <- average_cond(x75_imp, "Intensity.R_05_") / average_cond(x75_imp, "Intensity.C_05_")
fc_1_b1 <- average_cond(x75_imp, "Intensity.R_1_") / average_cond(x75_imp, "Intensity.C_1_")
fc_2_b1 <- average_cond(x75_imp, "Intensity.R_2_") / average_cond(x75_imp, "Intensity.C_2_")
fc_4_b1 <- average_cond(x75_imp, "Intensity.R_4_") / average_cond(x75_imp, "Intensity.C_4_")
fc_8_b1 <- average_cond(x75_imp, "Intensity.R_8_") / average_cond(x75_imp, "Intensity.C_8_")
fc_16_b1 <- average_cond(x75_imp, "Intensity.R_16_") / average_cond(x75_imp, "Intensity.C_16_")
fc_32_b1 <- average_cond(x75_imp, "Intensity.R_32_") / average_cond(x75_imp, "Intensity.C_32_")
fc_64_b1 <- average_cond(x75_imp, "Intensity.R_64_") / average_cond(x75_imp, "Intensity.C_64_")
fc_all_b1 <- cbind(fc_0_b1, fc_05_b1, fc_1_b1, fc_2_b1, fc_4_b1, fc_8_b1, fc_16_b1, fc_32_b1, fc_64_b1)

fc_0_b1_prot <- average_cond(x_prot50_imp, "Intensity[.]0_") / average_cond(x_prot50_imp, "Intensity[.]0_")
fc_05_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_05_") / average_cond(x_prot50_imp, "Intensity.C_05_")
fc_1_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_1_") / average_cond(x_prot50_imp, "Intensity.C_1_")
fc_2_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_2_") / average_cond(x_prot50_imp, "Intensity.C_2_")
fc_4_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_4_") / average_cond(x_prot50_imp, "Intensity.C_4_")
fc_8_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_8_") / average_cond(x_prot50_imp, "Intensity.C_8_")
fc_16_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_16_") / average_cond(x_prot50_imp, "Intensity.C_16_")
fc_32_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_32_") / average_cond(x_prot50_imp, "Intensity.C_32_")
fc_64_b1_prot <- average_cond(x_prot50_imp, "Intensity.R_64_") / average_cond(x_prot50_imp, "Intensity.C_64_")
fc_all_b1_prot <- cbind(fc_0_b1_prot, fc_05_b1_prot, fc_1_b1_prot, fc_2_b1_prot, fc_4_b1_prot, fc_8_b1_prot, fc_16_b1_prot, fc_32_b1_prot, fc_64_b1_prot)

########## Strategy 2 ##########
### 75 % quantile + median normalization (only for phosphoproteomics data)
# x <- x_old 
# x75 <- log2(as.matrix(x) + 1)
# 
# ind_batch1 <- grep(colnames(x75), pattern="0_1|05_1|1_1|2_1|4_1|8_1|16_1|32_1|64_1")
# ind_batch2 <- grep(colnames(x75), pattern="0_2|05_2|1_2|2_2|4_2|8_2|16_2|32_2|64_2")
# ind_batch3 <- grep(colnames(x75), pattern="0_3|05_3|1_3|2_3|4_3|8_3|16_3|32_3|64_3")
# ind_batch4 <- grep(colnames(x75), pattern="0_4|05_4|1_4|2_4|4_4|8_4|16_4|32_4|64_4")
# ind_batch5 <- grep(colnames(x75), pattern="0_5|05_5|1_5|2_5|4_5|8_5|16_5|32_5|64_5")
# 
# ind_batch_l <- list(ind_batch1, ind_batch2, ind_batch3) ##, ind_batch4, ind_batch5)
# 
# ## Day Normalization
# x_DNtemp <- x75
# for (i in 1:length(ind_batch_l)) {
#     subdata <- x[, ind_batch_l[[i]] ]
#     metab.med <- apply(subdata,1,median,na.rm=T)
#     x_DNtemp[, ind_batch_l[[i]] ] <- sweep(subdata,1,metab.med , FUN="-")
# }
# 
# ## add median per metabolite of whole experiment.
# x.med.all <- apply(x75 ,1, median, na.rm=T)
# x_DN <- x_DNtemp + x.med.all
# 
# ## impute
# ## only use replicates 1-3 for phospoproteomics
# x_DN <- x_DN[, -grep(colnames(x_DN), pattern=glob2rx("Intensity.?_*_4"))]
# x_DN <- x_DN[, -grep(colnames(x_DN), pattern=glob2rx("Intensity.?_*_5"))]
# x_DN <- x_DN[, -grep(colnames(x_DN), pattern=glob2rx("Intensity.0_4"))]
# x_DN <- x_DN[, -grep(colnames(x_DN), pattern=glob2rx("Intensity.0_5"))]
# 
# ## impute missing values (use min per row)
# x_DN_imp <- impute_min(x_DN)
# x_DN_imp <- shapiro_filter(x_DN_imp)
# 
# 
# ########## Strategy 3 ##########
# ## linear model to dect and remove batch effects (limma::removeBatchEffect)
# ##sample.med <- apply(x_DN,2, median, na.rm=T)
# ##sweep(x_DN,2, sample.med, FUN="-") + median(sample.med)
# 
# ## only use replicates 1-3 for phospoproteomics
# x_DMN <- x75[, -grep(colnames(x75), pattern=glob2rx("Intensity.?_*_4"))]
# x_DMN <- x_DMN[, -grep(colnames(x_DMN), pattern=glob2rx("Intensity.?_*_5"))]
# x_DMN <- x_DMN[, -grep(colnames(x_DMN), pattern=glob2rx("Intensity.0_4"))]
# x_DMN <- x_DMN[, -grep(colnames(x_DMN), pattern=glob2rx("Intensity.0_5"))]
# 
# 
# batch <- character(length(unlist(ind_batch_l)))
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.0_1"))] <- "red"
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.?_*_1"))] <- "red"
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.0_2"))] <- "blue"
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.?_*_2"))] <- "blue"
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.0_3"))] <- "green"
# batch[grep(colnames(x_DMN), pattern=glob2rx("Intensity.?_*_3"))] <- "green"
# ## apply function removeBatchEffect
# x_DMN <- limma::removeBatchEffect(x_DMN, batch=batch)
# 
# ## impute missing values (use min per row)
# x_DMN_imp <- impute_min(x_DMN)
# x_DMN_imp <- shapiro_filter(x_DMN_imp)
# 
# ################################
# 
# stopifnot(is.numeric(x_DN_imp))
# stopifnot(is.numeric(x_DMN_imp))
# 
# ################################################################################
# library(pheatmap)
# pheatmap(x_DN_imp, scale="row", show_rownames=FALSE)
# pheatmap(x_DMN_imp, scale="row", show_rownames=FALSE)
# 



# 
# rpcaDN <- pca(t(x_DN_imp), method="ppca", center=TRUE, scale="none", pPcs=5, seed=3455)
# slplot(rpcaDN, pcs=c(1,2), scoresLoadings=c(T, F), scol=batch)
# rpcaDMN <- pca(t(x_DMN_imp), method="ppca", center=TRUE, scale="none", pPcs=5, seed=3455)
# slplot(rpcaDMN, pcs=c(1,2), scoresLoadings=c(T, F), scol=batch)
# pheatmap(x_DMN_prot, scale="row") ## 4/5 cluster together 
################################################################################

## use DN for A2 and B2
#### A2 ####
## makes mathematically no sense to calculate fold change from batch median 
## normalized values (negative values)
# t_0_a2 <- x_DN[, grep(colnames(x_DN), pattern="Intensity.0_")]
# fc_0_a2 <- t_0_a2 / t_0_a2
# fc_05_a2 <- calculate_fc(x_DN, "^Intensity.*_05_")
# fc_1_a2 <- calculate_fc(x_DN, "^Intensity.*_1_") 
# fc_2_a2 <- calculate_fc(x_DN, "^Intensity.*_2_")
# fc_4_a2 <- calculate_fc(x_DN, "^Intensity.*_4_")
# fc_8_a2 <- calculate_fc(x_DN, "^Intensity.*_8_")
# fc_16_a2 <- calculate_fc(x_DN, "^Intensity.*_16_")
# fc_32_a2 <- calculate_fc(x_DN, "^Intensity.*_32_")
# fc_64_a2 <- calculate_fc(x_DN, "^Intensity.*_64_")

# fc_0_a2_m <- apply(fc_0_a2, 1, mean, na.rm=TRUE) 
# fc_05_a2_m <- apply(fc_05_a2, 1, mean, na.rm=TRUE)
# fc_1_a2_m <- apply(fc_1_a2, 1, mean, na.rm=TRUE)
# fc_2_a2_m <- apply(fc_2_a2, 1, mean, na.rm=TRUE)
# fc_4_a2_m <- apply(fc_4_a2, 1, mean, na.rm=TRUE)
# fc_8_a2_m <- apply(fc_8_a2, 1, mean, na.rm=TRUE)
# fc_16_a2_m <- apply(fc_16_a2, 1, mean, na.rm=TRUE)
# fc_32_a2_m <- apply(fc_32_a2, 1, mean, na.rm=TRUE)
# fc_64_a2_m <- apply(fc_64_a2, 1, mean, na.rm=TRUE)
# fc_all_a2 <- cbind(fc_0_a2_m, fc_05_a2_m, fc_1_a2_m, fc_2_a2_m, fc_4_a2_m, fc_8_a2_m, fc_16_a2_m, fc_32_a2_m, fc_64_a2_m)

#### B2 ####
## average first, then calculate fc
# fc_0_b2 <- average_cond(x_DN, "Intensity[.]0_") / average_cond(x_DN, "Intensity[.]0_")
# fc_05_b2 <- average_cond(x_DN, "Intensity.R_05_") / average_cond(x_DN, "Intensity.C_05_")
# fc_1_b2 <- average_cond(x_DN, "Intensity.R_1_") / average_cond(x_DN, "Intensity.C_1_")
# fc_2_b2 <- average_cond(x_DN, "Intensity.R_2_") / average_cond(x_DN, "Intensity.C_2_")
# fc_4_b2 <- average_cond(x_DN, "Intensity.R_4_") / average_cond(x_DN, "Intensity.C_4_")
# fc_8_b2 <- average_cond(x_DN, "Intensity.R_8_") / average_cond(x_DN, "Intensity.C_8_")
# fc_16_b2 <- average_cond(x_DN, "Intensity.R_16_") / average_cond(x_DN, "Intensity.C_16_")
# fc_32_b2 <- average_cond(x_DN, "Intensity.R_32_") / average_cond(x_DN, "Intensity.C_32_")
# fc_64_b2 <- average_cond(x_DN, "Intensity.R_64_") / average_cond(x_DN, "Intensity.C_64_")
# fc_all_b2 <- cbind(fc_0_b2, fc_05_b2, fc_1_b2, fc_2_b2, fc_4_b2, fc_8_b2, fc_16_b2, fc_32_b2, fc_64_b2)

## use DMN for A3 and B3
#### A3 ####
# t_0_a3 <- x_DMN[, grep(colnames(x_DMN), pattern="Intensity.0_")]
# fc_0_a3 <- t_0_a3 / t_0_a3
# fc_05_a3 <- calculate_fc(x_DMN, "^Intensity.*_05_")
# fc_1_a3 <- calculate_fc(x_DMN, "^Intensity.*_1_") 
# fc_2_a3 <- calculate_fc(x_DMN, "^Intensity.*_2_")
# fc_4_a3 <- calculate_fc(x_DMN, "^Intensity.*_4_")
# fc_8_a3 <- calculate_fc(x_DMN, "^Intensity.*_8_")
# fc_16_a3 <- calculate_fc(x_DMN, "^Intensity.*_16_")
# fc_32_a3 <- calculate_fc(x_DMN, "^Intensity.*_32_")
# fc_64_a3 <- calculate_fc(x_DMN, "^Intensity.*_64_")
# 
# fc_0_a3_m <- apply(fc_0_a3, 1, mean, na.rm=TRUE) 
# fc_05_a3_m <- apply(fc_05_a3, 1, mean, na.rm=TRUE)
# fc_1_a3_m <- apply(fc_1_a3, 1, mean, na.rm=TRUE)
# fc_2_a3_m <- apply(fc_2_a3, 1, mean, na.rm=TRUE)
# fc_4_a3_m <- apply(fc_4_a3, 1, mean, na.rm=TRUE)
# fc_8_a3_m <- apply(fc_8_a3, 1, mean, na.rm=TRUE)
# fc_16_a3_m <- apply(fc_16_a3, 1, mean, na.rm=TRUE)
# fc_32_a3_m <- apply(fc_32_a3, 1, mean, na.rm=TRUE)
# fc_64_a3_m <- apply(fc_64_a3, 1, mean, na.rm=TRUE)
# fc_all_a3 <- cbind(fc_0_a3_m, fc_05_a3_m, fc_1_a3_m, fc_2_a3_m, fc_4_a3_m, fc_8_a3_m, fc_16_a3_m, fc_32_a3_m, fc_64_a3_m)

#### B3 ####
## average first, then calculate fc
# fc_0_b3 <- average_cond(x_DMN, "Intensity[.]0_") / average_cond(x_DMN, "Intensity[.]0_")
# fc_05_b3 <- average_cond(x_DMN, "Intensity.R_05_") / average_cond(x_DMN, "Intensity.C_05_")
# fc_1_b3 <- average_cond(x_DMN, "Intensity.R_1_") / average_cond(x_DMN, "Intensity.C_1_")
# fc_2_b3 <- average_cond(x_DMN, "Intensity.R_2_") / average_cond(x_DMN, "Intensity.C_2_")
# fc_4_b3 <- average_cond(x_DMN, "Intensity.R_4_") / average_cond(x_DMN, "Intensity.C_4_")
# fc_8_b3 <- average_cond(x_DMN, "Intensity.R_8_") / average_cond(x_DMN, "Intensity.C_8_")
# fc_16_b3 <- average_cond(x_DMN, "Intensity.R_16_") / average_cond(x_DMN, "Intensity.C_16_")
# fc_32_b3 <- average_cond(x_DMN, "Intensity.R_32_") / average_cond(x_DMN, "Intensity.C_32_")
# fc_64_b3 <- average_cond(x_DMN, "Intensity.R_64_") / average_cond(x_DMN, "Intensity.C_64_")
# fc_all_b3 <- cbind(fc_0_b3, fc_05_b3, fc_1_b3, fc_2_b3, fc_4_b3, fc_8_b3, fc_16_b3, fc_32_b3, fc_64_b3)

##### visualisation ####
fc_all_a1_NAr <- fc_all_a1[apply(fc_all_a1, 1, function(x) !all(x == 1)),]
fc_all_b1_NAr <- fc_all_b1[apply(fc_all_b1, 1, function(x) !all(x == 1)),]
pheatmap(fc_all_a1_NAr, cluster_cols=FALSE, cluster_rows=TRUE, scale="row", show_rownames=F)
pheatmap(fc_all_b1_NAr, cluster_cols=FALSE, cluster_rows=TRUE, scale="row", show_rownames=F)

fc_all_a1_NAr_prot <- fc_all_a1_prot[apply(fc_all_a1_prot, 1, function(x) !all(x == 1)),]
fc_all_b1_NAr_prot <- fc_all_b1_prot[apply(fc_all_b1_prot, 1, function(x) !all(x == 1)),]
pheatmap(fc_all_a1_NAr_prot, cluster_cols=FALSE, cluster_rows=TRUE, scale="row", show_rownames=F)
pheatmap(fc_all_b1_NAr_prot, cluster_cols=FALSE, cluster_rows=TRUE, scale="row", show_rownames=F)

## create table with fold changes per sample per 
fc_x <- cbind(fc_0_a1, fc_05_a1, fc_1_a1, fc_2_a1, fc_4_a1, fc_8_a1, fc_16_a1, fc_32_a1, fc_64_a1)
fc_x_prot <- cbind(fc_0_a1_prot, fc_05_a1_prot, fc_1_a1_prot, fc_2_a1_prot, fc_4_a1_prot, fc_8_a1_prot, fc_16_a1_prot, fc_32_a1_prot, fc_64_a1_prot)

cols_x_sort <- c("Intensity.0_1", "Intensity.0_2", "Intensity.0_3", "Intensity.C_05_1",
    "Intensity.C_05_2", "Intensity.C_05_3", "Intensity.C_1_1", "Intensity.C_1_2",
    "Intensity.C_1_3",  "Intensity.C_2_1", "Intensity.C_2_2", "Intensity.C_2_3",
    "Intensity.C_4_1", "Intensity.C_4_2", "Intensity.C_4_3", "Intensity.C_8_1",
    "Intensity.C_8_2", "Intensity.C_8_3", "Intensity.C_16_1", "Intensity.C_16_2", 
    "Intensity.C_16_3", "Intensity.C_32_1", "Intensity.C_32_2", "Intensity.C_32_3",
    "Intensity.C_64_1", "Intensity.C_64_2", "Intensity.C_64_3", 
    "Intensity.0_1", "Intensity.0_2", "Intensity.0_3", "Intensity.R_05_1", 
    "Intensity.R_05_2", "Intensity.R_05_3", "Intensity.R_1_1", "Intensity.R_1_2",  
    "Intensity.R_1_3", "Intensity.R_2_1", "Intensity.R_2_2", "Intensity.R_2_3",
    "Intensity.R_4_1",  "Intensity.R_4_2",  "Intensity.R_4_3", "Intensity.R_8_1", 
    "Intensity.R_8_2", "Intensity.R_8_3", "Intensity.R_16_1", "Intensity.R_16_2", 
    "Intensity.R_16_3", "Intensity.R_32_1", "Intensity.R_32_2", "Intensity.R_32_3",
    "Intensity.R_64_1", "Intensity.R_64_2" ,"Intensity.R_64_3")

##x_DMN_sort <- x_DMN_imp[, cols_x_sort]
x_sort <- x75_imp[, cols_x_sort]
pheatmap(x_sort, scale="row", show_rownames=FALSE)
##pheatmap(x_DMN_sort, scale="row", show_rownames=FALSE)

x_prot_sort <- x_prot50_imp[, c("Intensity.0_1", "Intensity.0_2", "Intensity.0_3", "Intensity.0_4", "Intensity.0_5", 
    "Intensity.C_05_1", "Intensity.C_05_2", "Intensity.C_05_3",  "Intensity.C_05_4",  "Intensity.C_05_5", 
    "Intensity.C_1_1", "Intensity.C_1_2", "Intensity.C_1_3", "Intensity.C_1_4", "Intensity.C_1_5",  
    "Intensity.C_2_1", "Intensity.C_2_2", "Intensity.C_2_3", "Intensity.C_2_4", "Intensity.C_2_5",
    "Intensity.C_4_1", "Intensity.C_4_2", "Intensity.C_4_3", "Intensity.C_4_4", "Intensity.C_4_5", 
    "Intensity.C_8_1", "Intensity.C_8_2", "Intensity.C_8_3", "Intensity.C_8_4", "Intensity.C_8_5", 
    "Intensity.C_16_1", "Intensity.C_16_2", "Intensity.C_16_3", "Intensity.C_16_4", "Intensity.C_16_5", 
    "Intensity.C_32_1", "Intensity.C_32_2", "Intensity.C_32_3", "Intensity.C_32_4", "Intensity.C_32_5",
    "Intensity.C_64_1", "Intensity.C_64_2", "Intensity.C_64_3", "Intensity.C_64_4", "Intensity.C_64_5", 
    "Intensity.0_1", "Intensity.0_2", "Intensity.0_3", "Intensity.0_4", "Intensity.0_5", 
    "Intensity.R_05_1", "Intensity.R_05_2", "Intensity.R_05_3", "Intensity.R_05_4", "Intensity.R_05_5", 
    "Intensity.R_1_1", "Intensity.R_1_2", "Intensity.R_1_3", "Intensity.R_1_4", "Intensity.R_1_5", 
    "Intensity.R_2_1", "Intensity.R_2_2", "Intensity.R_2_3", "Intensity.R_2_4", "Intensity.R_2_5", 
    "Intensity.R_4_1",  "Intensity.R_4_2",  "Intensity.R_4_3", "Intensity.R_4_4", "Intensity.R_4_5", 
    "Intensity.R_8_1", "Intensity.R_8_2", "Intensity.R_8_3", "Intensity.R_8_4", "Intensity.R_8_5",
    "Intensity.R_16_1", "Intensity.R_16_2", "Intensity.R_16_3", "Intensity.R_16_4", "Intensity.R_16_5", 
    "Intensity.R_32_1", "Intensity.R_32_2", "Intensity.R_32_3", "Intensity.R_32_4", "Intensity.R_32_5",
    "Intensity.R_64_1", "Intensity.R_64_2" , "Intensity.R_64_3", "Intensity.R_64_4", "Intensity.R_64_5"
)]
pheatmap(x_prot_sort, scale="row")

## find significant features
## 2 factor anova
treatment <- c(rep("C", 3), rep(rep("C", 3), 8), rep("R", 3), rep(rep("R", 3), 8))
treatment_prot <- c(rep("C", 5), rep(rep("C", 5), 8), rep("R", 5), rep(rep("R", 5), 8))
time <- c(rep(c(rep(0, 3), rep(0.5, 3), rep(1, 3), rep(2,3), rep(4, 3), rep(8, 3), rep(16, 3), rep(32, 3), rep(64, 3)), 2))
time_fc <- c(rep(0, 3), rep(0.5, 3), rep(1, 3), rep(2, 3), rep(4, 3), rep(8, 3), rep(16, 3), rep(32, 3), rep(64, 3))
time_prot <- c(rep(c(rep(0, 5), rep(0.5, 5), rep(1, 5), rep(2, 5), rep(4, 5), rep(8, 5), rep(16, 5), rep(32, 5), rep(64, 5)), 2))
time_fc_prot <- c(rep(0, 5), rep(0.5, 5), rep(1, 5), rep(2, 5), rep(4, 5), rep(8, 5), rep(16, 5), rep(32, 5), rep(64, 5))

## overlap protein IDs of x_DMN (75% quantile + log2 + linear modell) and x_prot50 (50% quantile + log2)
## overlap protein IDs of x75 (75% quantile + log2) and x_prot50 (50% quantile + log2)
## continue analysis with three data sets 
## 1. x_prot50
## 2. x_DMN and x75
## 3. x_DMN/x_prot50 and x_DMN/x_prot50 that has corresponding ID --> normalize by protein content per replicate

## --> do the same for x with only 75% and log2 --> get overlap
## overlap: take ps --> significant, normalize per replicate prot
## whole proteome (not only !overlap): significant
## whole ps (not only !overlap): signicant

## check if rownames of x_DMN_sort and x_sort are identical
if (!all(rownames(x_DMN_sort) == rownames(x_sort))) stop()
id_phsi <- rownames(x_sort)
id_phsi <- unlist(lapply(strsplit(id_phsi, split="_"), "[", 1))
##id_phsi_DMN <- rownames(x_DMN_sort)
##id_phsi_DMN <- unlist(lapply(strsplit(id_phsi_DMN, split="_"), "[", 1))
id_prot <- rownames(x_prot_sort)

## 1. take all proteins
# inds_1 <- match(id_prot, id_phsi) 
# inds_1 <- is.na(inds_1)
x_prot_1 <- x_prot_sort ## [inds_1, ]
x_prot_fc_1 <- fc_x_prot ## [inds_1, ]

x_prot_fc_1 <- ifelse(x_prot_fc_1 == Inf, NA, x_prot_fc_1)
x_prot_fc_1 <- ifelse(x_prot_fc_1 == -Inf, NA, x_prot_fc_1)

## 2. take all phosphosites from x_DMN and x75
# inds_2 <- match(id_phsi, id_prot)
# inds_2 <- is.na(inds_2)
##x_DMN_2 <- x_DMN_sort # [inds_2, ]
x_2 <- x_sort
x_fc_2 <- fc_x

x_fc_2 <- ifelse(x_fc_2 == Inf, NA, x_fc_2)
x_fc_2 <- ifelse(x_fc_2 == -Inf, NA, x_fc_2)

## 3. 
inds_3 <- match(id_phsi, id_prot)
inds_3 <- !is.na(inds_3)
##inds_DMN_3 <- match(id_phsi_DMN, id_prot)
##inds_DMN_3 <- !is.na(inds_DMN_3)
x_3 <- x_sort[inds_3, ]
##x_DMN_3 <- x_DMN_sort[inds_DMN_3, ]
id_phsi_3 <- id_phsi[inds_3]
##id_phsi_DMN_3 <- id_phsi_DMN[inds_DMN_3]
## normalize phosphosite by protein content
inds_3_cut <- match(id_phsi_3, id_prot)
##inds_3_DMN_cut <- match(id_phsi_DMN_3, id_prot)
#x_DMN_norm_3 <- x_DMN_3 / x_prot_sort[inds_3_DMN_cut, colnames(x_prot_sort) %in% colnames(x_DMN_sort)]
x_norm_3 <- x_3 / x_prot_sort[inds_3_cut, colnames(x_prot_sort) %in% colnames(x_sort)]

x_norm_3 <- ifelse(x_norm_3 == Inf, NA, x_norm_3)
x_norm_3 <- ifelse(x_norm_3 == -Inf, NA, x_norm_3)


## calculate anova per feature and return error messages if not possible
aov_l <- function(x=x_DMN_sort, treatment=treatment, time=time) {
    aov_res <- lapply(seq_len(nrow(x)), function(y) {
        df <- data.frame(value=as.numeric(x[y,]), treatment=treatment, time=time)
        ## balanced design anova for treatment and time and interaction
        aov_res <- tryCatch(aov(value ~ treatment * time, data=df), error=function(i) print(y)) 
        return(aov_res)
    })
    return(aov_res)
}

aov_fc_l <- function(fc_x=fc_DMN, time_fc=time_fc) {
    aov_res <- lapply(seq_len(nrow(fc_x)), function(y) {
        df <- data.frame(value=as.numeric(fc_x[y,]), time=time_fc)
        ## balanced design anova for treatment and time and interaction
        aov_res <- tryCatch(aov(value ~ time, data=df), error=function(i) print(y))
        return(aov_res)
    })
    return(aov_res)
}

get_pvalue <- function(aov_l) {
    p_treatment <- lapply(seq_along(aov_l), function(x) try(summary(aov_l[[x]])[[1]]["Pr(>F)"][1, ]))
    p_time <- lapply(seq_along(aov_l), function(x) try(summary(aov_l[[x]])[[1]]["Pr(>F)"][2, ]))
    p_treatment_time <- lapply(seq_along(aov_l), function(x) try(summary(aov_l[[x]])[[1]]["Pr(>F)"][3, ]))
    ## unlist
    p_treatment <- as.numeric(unlist(p_treatment))
    p_time <- as.numeric(unlist(p_time))
    p_treatment_time <- as.numeric(unlist(p_treatment_time))
    
    return(list(p_treatment=p_treatment, p_time=p_time, 
                p_treatment_time=p_treatment_time))
}

get_pvalue_fc <- function(aov_fc_l) {
    p_time_fc <- lapply(seq_len(length(aov_fc_l)), function(x) try(summary(aov_fc_l[[x]])[[1]]["Pr(>F)"][1, ]))
    ## unlist
    p_time_fc <- as.numeric(unlist(p_time_fc))
    
    return(p_time_fc)
}

adj_pvalue <- function(p_value) {
    p_treatment_adj <- p.adjust(p_value[[1]], method="fdr")
    p_time_adj <- p.adjust(p_value[[2]], method="fdr")
    p_treatment_time_adj <- p.adjust(p_value[[3]], method="fdr")
    return(list(p_treatment_adj=p_treatment_adj, p_time_adj=p_time_adj, 
                p_treatment_time_adj=p_treatment_time_adj))
}

adj_pvalue_fc <- function(p_value_fc) {
    p_time_adj <- p.adjust(p_value_fc, method="fdr")
    return(p_time_adj)
}

## get p-values for treatment, time and interaction
## 1. 
aov_1 <- aov_l(x_prot_1, treatment_prot, time_prot)
aov_fc_1 <- aov_fc_l(x_prot_fc_1, time_fc_prot)
pvalue_1 <- get_pvalue(aov_1)
pvalue_adj_1 <- adj_pvalue(pvalue_1)
pvalue_fc_1 <- get_pvalue_fc(aov_fc_1)
pvalue_fc_adj_1 <- adj_pvalue_fc(pvalue_fc_1)

## 2. 
#aov_DMN_2 <- aov_l(x_DMN_2, treatment, time)
aov_2 <- aov_l(x_2, treatment, time)
aov_fc_2 <- aov_fc_l(x_fc_2, time_fc)
#pvalue_DMN_2 <- get_pvalue(aov_DMN_2)
pvalue_2 <- get_pvalue(aov_2)
pvalue_fc_2 <- get_pvalue_fc(aov_fc_2)
#pvalue_DMN_adj_2 <- adj_pvalue(pvalue_DMN_2)
pvalue_adj_2 <- adj_pvalue(pvalue_2)
pvalue_fc_adj_2 <- adj_pvalue_fc(pvalue_fc_2)

## 3.
#aov_DMN_norm_3 <- aov_l(x_DMN_norm_3, treatment, time)
aov_norm_3 <- aov_l(x_norm_3, treatment, time)
#pvalue_DMN_norm_3 <- get_pvalue(aov_DMN_norm_3)
pvalue_norm_3 <- get_pvalue(aov_norm_3)
#pvalue_DMN_norm_adj_3 <- adj_pvalue(pvalue_DMN_norm_3)
pvalue_norm_adj_3 <- adj_pvalue(pvalue_norm_3)

## fuzzy clustering and VS clustering
library(Mfuzz)
library(limma)
library(shiny)
library(qvalue)
library(matrixStats)
source("functions_clustering.R")

## 1. 
x_prot_1_sign <- x_prot_1[sort(unique(c(which(pvalue_adj_1[[1]] < 0.05), which(pvalue_adj_1[[2]] < 0.05), which(pvalue_adj_1[[3]] < 0.05)))), ]
x_prot_fc_1_sign <- x_prot_fc_1[rownames(x_prot_1_sign),]
## statWrapper
statOut_prot_fc_1 <- statWrapper(dat=x_prot_fc_1_sign, NumReps=5, NumCond=9, isPaired=F, isStat=T) ## only three feature
statOut_prot_1 <- statWrapper(dat=x_prot_1_sign, NumReps=5, NumCond=18, isPaired=F, isStat=T)
## estimClustNum
clustNumOut_prot_fc_1 <- estimClustNum(dat=statOut_prot_fc_1$dat, maxClust=20, cores=1)
clustNumOut_prot_1 <- estimClustNum(dat=statOut_prot_1$dat, maxClust=20, cores=1)
## runClustWrapper
ClustOut_vs_prot_fc_2 <- runClustWrapper(dat=statOut_prot_fc_1$dat, clustNumOut_prot_fc_1$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_prot_fc_2 <- runClustWrapper(dat=statOut_prot_fc_1$dat, clustNumOut_prot_fc_1$numclust_st, proteins=NULL, VSClust=F, cores=1)
ClustOut_vs_prot_1 <- runClustWrapper(dat=statOut_prot_1$dat, clustNumOut_prot_1$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_prot_1 <- runClustWrapper(dat=statOut_prot_1$dat, clustNumOut_prot_1$numclust_st, proteins=NULL, VSClust=F, cores=1)

## 2. 
x_2_sign <- x_2[sort(unique(c(which(pvalue_adj_2[[1]] < 0.05), which(pvalue_adj_2[[2]] < 0.05), which(pvalue_adj_2[[3]] < 0.05)))), ]
x_fc_2_sign <- x_fc_2[rownames(x_2_sign),]
##x_DMN_2_sign <- x_DMN_2[sort(unique(c(which(pvalue_DMN_adj_2[[1]] < 0.05), which(pvalue_DMN_adj_2[[2]] < 0.05), which(pvalue_DMN_adj_2[[3]] < 0.05)))), ]

## statWrapper
statOut_x_fc_2 <- statWrapper(dat=x_fc_2_sign, NumReps=3, NumCond=9, isPaired=F, isStat=T) ## ?? too few ??
##statOut_x_DMN_2 <- statWrapper(dat=x_DMN_2_sign, NumReps=3, NumCond=18, isPaired=F, isStat=T)
statOut_x_2 <- statWrapper(dat=x_2_sign, NumReps=3, NumCond=18, isPaired=F, isStat=T)
## estimClustNum
clustNumOut_x_fc_2 <- estimClustNum(dat=statOut_x_fc_2$dat, maxClust=20, cores=1)
##clustNumOut_x_DMN_2 <- estimClustNum(dat=statOut_x_DMN_2$dat, maxClust=20, cores=1)
clustNumOut_x_2 <- estimClustNum(dat=statOut_x_2$dat, maxClust=20, cores=1)
## run clustWrapper
ClustOut_vs_x_fc_2 <- runClustWrapper(dat=statOut_x_fc_2$dat, clustNumOut_x_fc_2$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_x_fc_2 <- runClustWrapper(dat=statOut_x_fc_2$dat, clustNumOut_x_fc_2$numclust_st, proteins=NULL, VSClust=F, cores=1)
#ClustOut_vs_x_DMN_2 <- runClustWrapper(dat=statOut_x_DMN_2$dat, clustNumOut_x_DMN_2$numclust_vs, proteins=NULL, VSClust=T, cores=1)
#ClustOut_st_x_DMN_2 <- runClustWrapper(dat=statOut_x_DMN_2$dat, clustNumOut_x_DMN_2$numclust_st, proteins=NULL, VSClust=F, cores=1)
ClustOut_vs_x_2 <- runClustWrapper(dat=statOut_x_2$dat, clustNumOut_x_2$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_x_2 <- runClustWrapper(dat=statOut_x_2$dat, clustNumOut_x_2$numclust_st, proteins=NULL, VSClust=F, cores=1)


## 3.
##x_DMN_norm_3_sign <- x_DMN_norm_3[sort(unique(c(which(pvalue_DMN_norm_adj_3[[1]] < 0.05), which(pvalue_DMN_norm_adj_3[[2]] < 0.05), which(pvalue_DMN_norm_adj_3[[3]] < 0.05)))),]
x_norm_3_sign <- x_norm_3[sort(unique(c(which(pvalue_norm_adj_3[[1]] < 0.05), which(pvalue_norm_adj_3[[2]] < 0.05), which(pvalue_norm_adj_3[[3]] < 0.05)))),]
## statWrapper
statOut_x_DMN_norm_3 <- statWrapper(dat=x_DMN_norm_3_sign, NumReps=3, NumCond=18, isPaired=F, isStat=T)
statOut_x_norm_3 <- statWrapper(dat=x_norm_3_sign, NumReps=3, NumCond=18, isPaired=F, isStat=T) ## ?? ## 
## estimClustNum
clustNumOut_x_DMN_norm_3 <- estimClustNum(dat=statOut_x_DMN_norm_3$dat, maxClust=20, cores=1)
clustNumOut_x_norm_3 <- estimClustNum(dat=statOut_x_norm_3$dat, maxClust=20, cores=1)
## runClustWrapper
ClustOut_vs_DMN_norm_3 <- runClustWrapper(dat=statOut_x_DMN_norm_3$dat, clustNumOut_x_DMN_norm_3$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_DMN_norm_3 <- runClustWrapper(dat=statOut_x_DMN_norm_3$dat, clustNumOut_x_DMN_norm_3$numclust_st, proteins=NULL, VSClust=F, cores=1)
ClustOut_vs_x_norm_3 <- runClustWrapper(dat=statOut_x_norm_3$dat, clustNumOut_x_norm_3$numclust_vs, proteins=NULL, VSClust=T, cores=1)
ClustOut_st_x_norm_3 <- runClustWrapper(dat=statOut_x_norm_3$dat, clustNumOut_x_norm_3$numclust_st, proteins=NULL, VSClust=F, cores=1)

if (normalization == "batchEffect") {
    write.table(x_prot_1_sign, file="x_prot_1_sign_batchEffect.txt", sep="\t", dec=".")
    write.table(x_2_sign, file="x_2_sign_batchEffect.txt", sep="\t", dec=".")
    write.table(x_norm_3_sign, file="x_norm_3_sign_batchEffect.txt", sep="\t", dec=".")
    write.table(x_prot_fc_1_sign, file="x_prot_fc_1_sign_batchEffect.txt", sep="\t", dec=".")
    write.table(x_fc_2_sign, file="x_fc_2_sign_batchEffect.txt", sep="\t", dec=".")
    
}

if (normalization == "quantile") {
    write.table(x_prot_1_sign, file="x_prot_1_sign_quantile.txt", sep="\t", dec=".")
    write.table(x_2_sign, file="x_2_sign_quantile.txt", sep="\t", dec=".")
    write.table(x_norm_3_sign, file="x_norm_3_sign_quantile.txt", sep="\t", dec=".")
    write.table(x_prot_fc_1_sign, file="x_prot_fc_1_sign_quantile.txt", sep="\t", dec=".")
    write.table(x_fc_2_sign, file="x_fc_2_sign_quantile.txt", sep="\t", dec=".")
}
