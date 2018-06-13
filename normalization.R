###### Phosphosites #####

setwd("H:/Documents/05_collaborations_small_projects/Umarah_phosphoproteomics_normalization")
x <- read.table("valid_filt_PhosphoProteomics_0_64_copy.txt", header = TRUE, fill = TRUE, sep = "\t")
x <- x[,grep(colnames(x), pattern = "Intensity")]

x <- apply(x, 2, function(x) x / quantile(x, .75, na.rm=T))
x <- x[, !colnames(x) %in% c("Intensity",  "Intensity.0_6", "Intensity.0_7")]

## check for outliers --> PCA and manually remove
library(pcaMethods)
x_outl <- x
x_outl <- ifelse(is.na(x_outl), NA, x)
rpcaNN <- pca(t(x_outl), method="ppca", center=TRUE, scale="none", nPcs=5, seed = 3455, na.rm = TRUE)
slplot(rpcaNN, pcs=c(1,2), scoresLoadings = c(T, F), scex=0.5) ## looks ok, do not remove outliers

## filtering by timepoints, if 3/5 we take it
inds <- grep(colnames(x), pattern = "Intensity[.]0_") ## 1|2, 3|4, 5|6, 7, 8
t0_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
t0_2 <- !apply(x[,inds[3:4]], 1, function(x) all(is.na(x)))
t0_3 <- !apply(x[,inds[5:6]], 1, function(x) all(is.na(x)))
t0_4 <- !is.na(x[,inds[7]])
t0_5 <- !is.na(x[,inds[8]]) 
t0 <- t0_1 + t0_2 + t0_3 + t0_4 + t0_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_0_5_") ## 1|2, 3|4, 5, 6, 7
tC05_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tC05_2 <- !apply(x[,inds[3:4]], 1, function(x) all(is.na(x)))
tC05_3 <- !is.na(x[,inds[5]])
tC05_4 <- !is.na(x[,inds[6]])
tC05_5 <- !is.na(x[,inds[7]]) 
tC05 <- tC05_1 + tC05_2 + tC05_3 + tC05_4 + tC05_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_0_5_") ## 1|2, 3, 4, 5, 6
tR05_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tR05_2 <- !is.na(x[,inds[3]])
tR05_3 <- !is.na(x[,inds[4]])
tR05_4 <- !is.na(x[,inds[5]])
tR05_5 <- !is.na(x[,inds[6]]) 
tR05 <- tR05_1 + tR05_2 + tR05_3 + tR05_4 + tR05_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_1_") ## 1|2, 3, 4, 5, 6
tC1_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tC1_2 <- !is.na(x[,inds[3]])
tC1_3 <- !is.na(x[,inds[4]])
tC1_4 <- !is.na(x[,inds[5]])
tC1_5 <- !is.na(x[,inds[6]]) 
tC1 <- tC1_1 + tC1_2 + tC1_3 + tC1_4 + tC1_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_1_") ## 1|2, 3, 4, 5, 6
tR1_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tR1_2 <- !is.na(x[,inds[3]])
tR1_3 <- !is.na(x[,inds[4]])
tR1_4 <- !is.na(x[,inds[5]])
tR1_5 <- !is.na(x[,inds[6]]) 
tR1 <- tR1_1 + tR1_2 + tR1_3 + tR1_4 + tR1_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_2_") ## 1|2, 3|4, 5, 6, 7
tC2_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tC2_2 <- !apply(x[,inds[3:4]], 1, function(x) all(is.na(x)))
tC2_3 <- !is.na(x[,inds[5]])
tC2_4 <- !is.na(x[,inds[6]])
tC2_5 <- !is.na(x[,inds[7]]) 
tC2 <- tC2_1 + tC2_2 + tC2_3 + tC2_4 + tC2_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_2_") ## 1|2, 3|4, 4, 5, 6
tR2_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tR2_2 <- !apply(x[,inds[3:4]], 1, function(x) all(is.na(x)))
tR2_3 <- !is.na(x[,inds[5]])
tR2_4 <- !is.na(x[,inds[6]])
tR2_5 <- !is.na(x[,inds[7]]) 
tR2 <- tR2_1 + tR2_2 + tR2_3 + tR2_4 + tR2_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_4_") ## 1|2, 3, 4, 5, 6
tC4_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tC4_2 <- !is.na(x[,inds[3]])
tC4_3 <- !is.na(x[,inds[4]])
tC4_4 <- !is.na(x[,inds[5]])
tC4_5 <- !is.na(x[,inds[6]]) 
tC4 <- tC4_1 + tC4_2 + tC4_3 + tC4_4 + tC4_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_4_") ## 1|2, 3, 4, 5, 6
tR4_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tR4_2 <- !is.na(x[,inds[3]])
tR4_3 <- !is.na(x[,inds[4]])
tR4_4 <- !is.na(x[,inds[5]])
tR4_5 <- !is.na(x[,inds[6]]) 
tR4 <- tR4_1 + tR4_2 + tR4_3 + tR4_4 + tR4_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_8_") ## 1, 2, 3, 4, 5
tC8_1 <- !is.na(x[,inds[1]])
tC8_2 <- !is.na(x[,inds[2]])
tC8_3 <- !is.na(x[,inds[3]])
tC8_4 <- !is.na(x[,inds[4]])
tC8_5 <- !is.na(x[,inds[5]]) 
tC8 <- tC8_1 + tC8_2 + tC8_3 + tC8_4 + tC8_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_8_") ## 1|2, 3, 4, NA, 5
tR8_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tR8_2 <- !is.na(x[,inds[3]])
tR8_3 <- !is.na(x[,inds[4]])
##tR8_4 <- !is.na(x[,inds[5]])
tR8_5 <- !is.na(x[,inds[5]]) 
tR8 <- tR8_1 + tR8_2 + tR8_3 + tR8_5 >= 2
inds <- grep(colnames(x), pattern = "Intensity[.]C_16_") ## 1|2|3, 4, 5, 6, 7
tC16_1 <- !apply(x[,inds[1:3]], 1, function(x) all(is.na(x)))
tC16_2 <- !is.na(x[,inds[4]])
tC16_3 <- !is.na(x[,inds[5]])
tC16_4 <- !is.na(x[,inds[6]])
tC16_5 <- !is.na(x[,inds[7]]) 
tC16 <- tC16_1 + tC16_2 + tC16_3 + tC16_4 + tC16_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_16_") ## 1, 2|3, 4, 5, 6
tR16_1 <- !is.na(x[,inds[1]])
tR16_2 <- !apply(x[,inds[2:3]], 1, function(x) all(is.na(x)))
tR16_3 <- !is.na(x[,inds[4]])
tR16_4 <- !is.na(x[,inds[5]])
tR16_5 <- !is.na(x[,inds[6]]) 
tR16 <- tR16_1 + tR16_2 + tR16_3 + tR16_4 + tR16_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_32_") ## 1|2, 3, 4, 5|6, NA
tC32_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x)))
tC32_2 <- !is.na(x[,inds[3]])
tC32_3 <- !is.na(x[,inds[4]])
tC32_4 <- !apply(x[,inds[5:6]], 1, function(x) all(is.na(x)))
##tC16_5 <- !is.na(x[,inds[]]) 
tC32 <- tC32_1 + tC32_2 + tC32_3 + tC32_4 >= 2
inds <- grep(colnames(x), pattern = "Intensity[.]R_32_") ## 1, 2, 3, 4, 5|6
tR32_1 <- !is.na(x[,inds[1]])
tR32_2 <- !is.na(x[,inds[2]])
tR32_3 <- !is.na(x[,inds[3]])
tR32_4 <- !is.na(x[,inds[4]])
tR32_5 <- !apply(x[,inds[5:6]], 1, function(x) all(is.na(x))) 
tR32 <- tR32_1 + tR32_2 + tR32_3 + tR32_4 + tR32_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]C_64_") ## 1, 2, 3, 4, 5
tC64_1 <- !is.na(x[,inds[1]])
tC64_2 <- !is.na(x[,inds[2]])
tC64_3 <- !is.na(x[,inds[3]])
tC64_4 <- !is.na(x[,inds[4]])
tC64_5 <- !is.na(x[,inds[5]]) 
tC64 <- tC64_1 + tC64_2 + tC64_3 + tC64_4 + tC64_5 >= 3
inds <- grep(colnames(x), pattern = "Intensity[.]R_64_") ## 1|2, 3|4, 5, 6, 7
tR64_1 <- !apply(x[,inds[1:2]], 1, function(x) all(is.na(x))) 
tR64_2 <- !apply(x[,inds[3:4]], 1, function(x) all(is.na(x))) 
tR64_3 <- !is.na(x[,inds[5]])
tR64_4 <- !is.na(x[,inds[6]])
tR64_5 <- !is.na(x[,inds[7]])
tR64 <- tR64_1 + tR64_2 + tR64_3 + tR64_4 + tR64_5 >= 3

## take these lines that are present in at least one of the timepoints with >= 3
t_keep <- t0 + tC05 + tR05 + tC1 + tR1 + tC2 + tR2 + tC4 + tR4 + tC8 + tR8 + tC16 + tR16 + tC32 + tR32 + tC64 + tR64 > 0
x <- x[t_keep, ]
x_outl <- x
x_outl <- ifelse(is.na(x_outl), NA, x)
rpcaNN <- pca(t(x_outl), method="ppca", center=TRUE, scale="none", nPcs=5, seed = 3455, na.rm = TRUE)
slplot(rpcaNN, pcs=c(1,2), scoresLoadings = c(T, F), scex=0.5) ## looks ok, do not remove outliers

## calculating average of replicates and delete unnecessary replicates
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_1")]; x <- cbind(x, "Intensity.0_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_2")]; x <- cbind(x, "Intensity.0_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_3")]; x <- cbind(x, "Intensity.0_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_1")]; x <- cbind(x, "Intensity.C_0_5_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_2")]; x <- cbind(x, "Intensity.C_0_5_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_1_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_1_1")]; x <- cbind(x, "Intensity.C_1_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_16_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_16_1")]; x <- cbind(x, "Intensity.C_16_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_1")]; x <- cbind(x, "Intensity.C_2_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_2")]; x <- cbind(x, "Intensity.C_2_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_32_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_32_1")]; x <- cbind(x, "Intensity.C_32_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_32_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_32_4")]; x <- cbind(x, "Intensity.C_32_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_4_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_4_1")]; x <- cbind(x, "Intensity.C_4_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_0_5_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_0_5_1")]; x <- cbind(x, "Intensity.R_0_5_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_1_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_1_1")]; x <- cbind(x, "Intensity.R_1_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_16_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_16_2")]; x <- cbind(x, "Intensity.R_16_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_2_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_2_1")]; x <- cbind(x, "Intensity.R_2_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_2_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_2_2")]; x <- cbind(x, "Intensity.R_2_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_32_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_32_5")]; x <- cbind(x, "Intensity.R_32_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_4_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_4_1")]; x <- cbind(x, "Intensity.R_4_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_64_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_64_1")]; x <- cbind(x, "Intensity.R_64_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_64_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_64_2")]; x <- cbind(x, "Intensity.R_64_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_8_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_8_1")]; x <- cbind(x, "Intensity.R_8_1" = tmp)
colnames(x)[which(colnames(x) == "Intensity.C_4_2_180225173818")] <-  "Intensity.C_4_2"
colnames(x) <- gsub(colnames(x), pattern = "C_0_5_", replacement = "C_05_") ## rename columns
colnames(x) <- gsub(colnames(x), pattern = "R_0_5_", replacement = "R_05_") ## rename columns

x <- x[, !colnames(x) %in% "Intensity.C_8_4"] ## remove since there is no Intensity.R_8_4
x <- x[, !colnames(x) %in% "Intensity.R_32_5"] ## remove since there is no Intensity.C_32_5
x <- x[, sort(colnames(x))]

## save x in x_old
x_old <- x 


## define functions that are used later
## calculate_fc is used in the strategies A
calculate_fc <- function(x, pattern = "^Intensity.*_05_") {
    ind <- grep(colnames(x), pattern = pattern)
    beg <- seq(1, length(ind) / 2) 
    end <- beg + length(ind) / 2
    fc <- x[, ind[end]] / x[, ind[beg]]
    return(fc)
} 
## average_cond is used in the strategies B
average_cond <- function(x, condition) {
    x_cond <- x[, grep(colnames(x), pattern = condition)]
    x_cond_mean <- apply(x_cond, 1, mean, na.rm = TRUE)
    return(x_cond_mean)
}

## Strategy 1 (only 75% quantile)
######### 1A  ##########
## fc per replicate
x <- x_old 
x <- log2(x+1)

values0 <- x[, grep(colnames(x), pattern = "Intensity.0_")]
fc_0_a1 <- values0 / values0
fc_05_a1 <- calculate_fc(x, "^Intensity.*_05_")
fc_1_a1 <- calculate_fc(x, "^Intensity.*_1_") 
fc_2_a1 <- calculate_fc(x, "^Intensity.*_2_")
fc_4_a1 <- calculate_fc(x, "^Intensity.*_4_")
fc_8_a1 <- calculate_fc(x, "^Intensity.*_8_")
fc_16_a1 <- calculate_fc(x, "^Intensity.*_16_")
fc_32_a1 <- calculate_fc(x, "^Intensity.*_32_")
fc_64_a1 <- calculate_fc(x, "^Intensity.*_64_")

fc_0_a1_m <- apply(fc_0_a1, 1, mean, na.rm = TRUE)
fc_05_a1_m <- apply(fc_05_a1, 1, mean, na.rm = TRUE)
fc_1_a1_m <- apply(fc_1_a1, 1, mean, na.rm = TRUE)
fc_2_a1_m <- apply(fc_2_a1, 1, mean, na.rm = TRUE)
fc_4_a1_m <- apply(fc_4_a1, 1, mean, na.rm = TRUE)
fc_8_a1_m <- apply(fc_8_a1, 1, mean, na.rm = TRUE)
fc_16_a1_m <- apply(fc_16_a1, 1, mean, na.rm = TRUE)
fc_32_a1_m <- apply(fc_32_a1, 1, mean, na.rm = TRUE)
fc_64_a1_m <- apply(fc_64_a1, 1, mean, na.rm = TRUE)
fc_all_a1 <- cbind(fc_0_a1_m, fc_05_a1_m, fc_1_a1_m, fc_2_a1_m, fc_4_a1_m, fc_8_a1_m, fc_16_a1_m, fc_32_a1_m, fc_64_a1_m)

##### Strategy 1 B-1 average first replicates then calculate fc ######
x <- x_old
x <- log2(x+1)

fc_0_b1 <- average_cond(x, "Intensity[.]0_") / average_cond(x_old, "Intensity[.]0_")
fc_05_b1 <- average_cond(x, "Intensity.R_05_") / average_cond(x_old, "Intensity.C_05_")
fc_1_b1 <- average_cond(x, "Intensity.R_1_") / average_cond(x_old, "Intensity.C_1_")
fc_2_b1 <- average_cond(x, "Intensity.R_2_") / average_cond(x_old, "Intensity.C_2_")
fc_4_b1 <- average_cond(x, "Intensity.R_4_") / average_cond(x_old, "Intensity.C_4_")
fc_8_b1 <- average_cond(x, "Intensity.R_8_") / average_cond(x_old, "Intensity.C_8_")
fc_16_b1 <- average_cond(x, "Intensity.R_16_") / average_cond(x_old, "Intensity.C_16_")
fc_32_b1 <- average_cond(x, "Intensity.R_32_") / average_cond(x_old, "Intensity.C_32_")
fc_64_b1 <- average_cond(x, "Intensity.R_64_") / average_cond(x_old, "Intensity.C_64_")
fc_all_b1 <- cbind(fc_0_b1, fc_05_b1, fc_1_b1, fc_2_b1, fc_4_b1, fc_8_b1, fc_16_b1, fc_32_b1, fc_64_b1)

########## Strategy 2 & 3 ##########
### 75 % quantile + median normalization 
x <- x_old 
ind_batch1 <- grep(colnames(x), pattern = "0_1|05_1|1_1|2_1|4_1|8_1|16_1|32_1|64_1")
ind_batch2 <- grep(colnames(x), pattern = "0_2|05_2|1_2|2_2|4_2|8_2|16_2|32_2|64_2")
ind_batch3 <- grep(colnames(x), pattern = "0_3|05_3|1_3|2_3|4_3|8_3|16_3|32_3|64_3")
ind_batch4 <- grep(colnames(x), pattern = "0_4|05_4|1_4|2_4|4_4|8_4|16_4|32_4|64_4")
ind_batch5 <- grep(colnames(x), pattern = "0_5|05_5|1_5|2_5|4_5|8_5|16_5|32_5|64_5")

## remove the rows that have more than 4 missing values per row
# ind_remove <- as.numeric(which(apply(x[, ind_batch1], 1, function(x) sum(is.na(x))) > 5))
# ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch2], 1, function(x) sum(is.na(x))) > 5)))
# ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch3], 1, function(x) sum(is.na(x))) > 5)))
# ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch4], 1, function(x) sum(is.na(x))) > 5)))
# ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch5], 1, function(x) sum(is.na(x))) > 5)))
# ind_remove <- unique(ind_remove)
# x <- x[-ind_remove, ]

ind_batch_l <- list(ind_batch1, ind_batch2, ind_batch3, ind_batch4, ind_batch5)

## Day Normalization
x_DNtemp <- x
for (i in 1:length(ind_batch_l)) {
    subdata <- x[, ind_batch_l[[i]] ]
    metab.med <- apply(subdata,1,median,na.rm=T)
    x_DNtemp[, ind_batch_l[[i]] ] <- sweep(subdata,1,metab.med , FUN="-")
}
## add median per metabolite of whole experiment.
x.med.all <- apply(x ,1, median, na.rm=T)
x_DN <- x_DNtemp + x.med.all

################################################################################
## sample median normalization
sample.med <- apply(x_DN,2, median, na.rm=T)
x_DMN <- sweep(x_DN,2, sample.med, FUN="-") + median(sample.med)

x_DN <- log2(as.matrix(x_DN) + 1)
x_DMN <- log2(as.matrix(x_DMN) + 1)

stopifnot(is.numeric(x_DN))
stopifnot(is.numeric(x_DMN))

x_DN <- ifelse(is.na(x_DN), NA, x_DN)
x_DMN <- ifelse(is.na(x_DMN), NA, x_DMN)

## use DN for A2 and B2
#### A2 ####
t_0_a2 <- x_DN[, grep(colnames(x_DN), pattern = "Intensity.0_")]
fc_0_a2 <- t_0_a2 / t_0_a2
fc_05_a2 <- calculate_fc(x_DN, "^Intensity.*_05_")
fc_1_a2 <- calculate_fc(x_DN, "^Intensity.*_1_") 
fc_2_a2 <- calculate_fc(x_DN, "^Intensity.*_2_")
fc_4_a2 <- calculate_fc(x_DN, "^Intensity.*_4_")
fc_8_a2 <- calculate_fc(x_DN, "^Intensity.*_8_")
fc_16_a2 <- calculate_fc(x_DN, "^Intensity.*_16_")
fc_32_a2 <- calculate_fc(x_DN, "^Intensity.*_32_")
fc_64_a2 <- calculate_fc(x_DN, "^Intensity.*_64_")

fc_0_a2_m <- apply(fc_0_a2, 1, mean, na.rm = TRUE) 
fc_05_a2_m <- apply(fc_05_a2, 1, mean, na.rm = TRUE)
fc_1_a2_m <- apply(fc_1_a2, 1, mean, na.rm = TRUE)
fc_2_a2_m <- apply(fc_2_a2, 1, mean, na.rm = TRUE)
fc_4_a2_m <- apply(fc_4_a2, 1, mean, na.rm = TRUE)
fc_8_a2_m <- apply(fc_8_a2, 1, mean, na.rm = TRUE)
fc_16_a2_m <- apply(fc_16_a2, 1, mean, na.rm = TRUE)
fc_32_a2_m <- apply(fc_32_a2, 1, mean, na.rm = TRUE)
fc_64_a2_m <- apply(fc_64_a2, 1, mean, na.rm = TRUE)
fc_all_a2 <- cbind(fc_0_a2_m, fc_05_a2_m, fc_1_a2_m, fc_2_a2_m, fc_4_a2_m, fc_8_a2_m, fc_16_a2_m, fc_32_a2_m, fc_64_a2_m)

#### B2 ####
## average first, then calculate fc
fc_0_b2 <- average_cond(x_DN, "Intensity[.]0_") / average_cond(x_DN, "Intensity[.]0_")
fc_05_b2 <- average_cond(x_DN, "Intensity.R_05_") / average_cond(x_DN, "Intensity.C_05_")
fc_1_b2 <- average_cond(x_DN, "Intensity.R_1_") / average_cond(x_DN, "Intensity.C_1_")
fc_2_b2 <- average_cond(x_DN, "Intensity.R_2_") / average_cond(x_DN, "Intensity.C_2_")
fc_4_b2 <- average_cond(x_DN, "Intensity.R_4_") / average_cond(x_DN, "Intensity.C_4_")
fc_8_b2 <- average_cond(x_DN, "Intensity.R_8_") / average_cond(x_DN, "Intensity.C_8_")
fc_16_b2 <- average_cond(x_DN, "Intensity.R_16_") / average_cond(x_DN, "Intensity.C_16_")
fc_32_b2 <- average_cond(x_DN, "Intensity.R_32_") / average_cond(x_DN, "Intensity.C_32_")
fc_64_b2 <- average_cond(x_DN, "Intensity.R_64_") / average_cond(x_DN, "Intensity.C_64_")
fc_all_b2 <- cbind(fc_0_b2, fc_05_b2, fc_1_b2, fc_2_b2, fc_4_b2, fc_8_b2, fc_16_b2, fc_32_b2, fc_64_b2)

## use DMN for A3 and B3
#### A3 ####
t_0_a3 <- x_DMN[, grep(colnames(x_DMN), pattern = "Intensity.0_")]
fc_0_a3 <- t_0_a3 / t_0_a3
fc_05_a3 <- calculate_fc(x_DMN, "^Intensity.*_05_")
fc_1_a3 <- calculate_fc(x_DMN, "^Intensity.*_1_") 
fc_2_a3 <- calculate_fc(x_DMN, "^Intensity.*_2_")
fc_4_a3 <- calculate_fc(x_DMN, "^Intensity.*_4_")
fc_8_a3 <- calculate_fc(x_DMN, "^Intensity.*_8_")
fc_16_a3 <- calculate_fc(x_DMN, "^Intensity.*_16_")
fc_32_a3 <- calculate_fc(x_DMN, "^Intensity.*_32_")
fc_64_a3 <- calculate_fc(x_DMN, "^Intensity.*_64_")

fc_0_a3_m <- apply(fc_0_a3, 1, mean, na.rm = TRUE) 
fc_05_a3_m <- apply(fc_05_a3, 1, mean, na.rm = TRUE)
fc_1_a3_m <- apply(fc_1_a3, 1, mean, na.rm = TRUE)
fc_2_a3_m <- apply(fc_2_a3, 1, mean, na.rm = TRUE)
fc_4_a3_m <- apply(fc_4_a3, 1, mean, na.rm = TRUE)
fc_8_a3_m <- apply(fc_8_a3, 1, mean, na.rm = TRUE)
fc_16_a3_m <- apply(fc_16_a3, 1, mean, na.rm = TRUE)
fc_32_a3_m <- apply(fc_32_a3, 1, mean, na.rm = TRUE)
fc_64_a3_m <- apply(fc_64_a3, 1, mean, na.rm = TRUE)
fc_all_a3 <- cbind(fc_0_a3_m, fc_05_a3_m, fc_1_a3_m, fc_2_a3_m, fc_4_a3_m, fc_8_a3_m, fc_16_a3_m, fc_32_a3_m, fc_64_a3_m)

#### B3 ####
## average first, then calculate fc
fc_0_b3 <- average_cond(x_DMN, "Intensity[.]0_") / average_cond(x_DMN, "Intensity[.]0_")
fc_05_b3 <- average_cond(x_DMN, "Intensity.R_05_") / average_cond(x_DMN, "Intensity.C_05_")
fc_1_b3 <- average_cond(x_DMN, "Intensity.R_1_") / average_cond(x_DMN, "Intensity.C_1_")
fc_2_b3 <- average_cond(x_DMN, "Intensity.R_2_") / average_cond(x_DMN, "Intensity.C_2_")
fc_4_b3 <- average_cond(x_DMN, "Intensity.R_4_") / average_cond(x_DMN, "Intensity.C_4_")
fc_8_b3 <- average_cond(x_DMN, "Intensity.R_8_") / average_cond(x_DMN, "Intensity.C_8_")
fc_16_b3 <- average_cond(x_DMN, "Intensity.R_16_") / average_cond(x_DMN, "Intensity.C_16_")
fc_32_b3 <- average_cond(x_DMN, "Intensity.R_32_") / average_cond(x_DMN, "Intensity.C_32_")
fc_64_b3 <- average_cond(x_DMN, "Intensity.R_64_") / average_cond(x_DMN, "Intensity.C_64_")
fc_all_b3 <- cbind(fc_0_b3, fc_05_b3, fc_1_b3, fc_2_b3, fc_4_b3, fc_8_b3, fc_16_b3, fc_32_b3, fc_64_b3)


##### visualisation ####
library(pheatmap)

fc_all_a1_NAr <- fc_all_a1[apply(fc_all_a1, 1, function(x) sum(is.na(x))) < 5,]
fc_all_b1_NAr <- fc_all_b1[apply(fc_all_b1, 1, function(x) sum(is.na(x))) < 5,]
fc_all_a2_NAr <- fc_all_a2[apply(fc_all_a2, 1, function(x) sum(is.na(x))) < 5,]
fc_all_b2_NAr <- fc_all_b2[apply(fc_all_b2, 1, function(x) sum(is.na(x))) < 5,]
fc_all_a3_NAr <- fc_all_a3[apply(fc_all_a3, 1, function(x) sum(is.na(x))) < 5,]
fc_all_b3_NAr <- fc_all_b3[apply(fc_all_b3, 1, function(x) sum(is.na(x))) < 5,]
pheatmap(fc_all_a1_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap(fc_all_b1_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap(fc_all_a2_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap(fc_all_b2_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap(fc_all_a3_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
pheatmap(fc_all_b3_NAr, cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")

# 
# rpcaNN <- pca(t(x), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455, na.rm = TRUE)
# rpcaDN <- pca(t(x_DN), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)
# rpcaDMN <- pca(t(x_DMN), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)
# rpcaQua <- pca(t(x_qu), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)
# 
# # define colors (experiment dependent)
# cls <- unlist(lapply(strsplit(colnames(x), split = "_"), function(x) paste(x[1], x[2], sep = "_")) )
# cls[grep(cls, pattern = "Intensity[.]0_")] <- "Intensity_0"
# cls <- unlist(lapply(strsplit(cls, split = "_"), "[", 2))
# cls <- as.factor(cls)
# cls <- as.numeric(cls)
# 
# #pchCR
# pchCR <- numeric(length(cls))
# pchCR[grep(colnames(x), pattern = "Intensity[.]C")] <- 1
# pchCR[grep(colnames(x), pattern = "Intensity[.]R")] <- 2
# pchCR[grep(colnames(x), pattern = "Intensity[.]0")] <- 3
# 
# #pdf('PCA_Mix_sample_NN.pdf',7,7)
# slplot(rpcaNN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
# slplot(rpcaNN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
# slplot(rpcaNN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
# #dev.off()
# 
# #pdf('PCA_Mix_sample_DN.pdf',7,7)
# slplot(rpcaDN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=0.5, spch = pchCR )
# slplot(rpcaDN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
# slplot(rpcaDN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
# #dev.off()
# 
# 
# #pdf('PCA_Mix_sample_DMN.pdf',7,7)
# slplot(rpcaDMN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
# slplot(rpcaDMN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
# slplot(rpcaDMN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
# #dev.off()
# 
# slplot(rpcaQua, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )







####################################################################################################
setwd("H:/Documents/05_collaborations_small_projects/Umarah_phosphoproteomics_normalization")
x <- read.table("valid_filt_Proteomics_0_64_copy.txt", header = TRUE, fill = TRUE, sep = "\t")
x <- x[,grep(colnames(x), pattern = "Intensity")]
## calculating average and delete
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_2")]; x <- cbind(x, "Intensity.0_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_3")]; x <- cbind(x, "Intensity.0_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_4")]; x <- cbind(x, "Intensity.0_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_6")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_6")]; x <- cbind(x, "Intensity.0_6" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]0_8")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]0_8")]; x <- cbind(x, "Intensity.0_8" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_2")]; x <- cbind(x, "Intensity.C_0_5_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_3")]; x <- cbind(x, "Intensity.C_0_5_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_4")]; x <- cbind(x, "Intensity.C_0_5_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_0_5_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_0_5_5")]; x <- cbind(x, "Intensity.C_0_5_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_1_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_1_2")]; x <- cbind(x, "Intensity.C_1_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_1_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_1_4")]; x <- cbind(x, "Intensity.C_1_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_1_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_1_5")]; x <- cbind(x, "Intensity.C_1_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_16_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_16_1")]; x <- cbind(x, "Intensity.C_16_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_16_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_16_2")]; x <- cbind(x, "Intensity.C_16_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_16_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_16_3")]; x <- cbind(x, "Intensity.C_16_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_1")]; x <- cbind(x, "Intensity.C_2_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_2")]; x <- cbind(x, "Intensity.C_2_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_3")]; x <- cbind(x, "Intensity.C_2_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_2_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_2_4")]; x <- cbind(x, "Intensity.C_2_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_32_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_32_1")]; x <- cbind(x, "Intensity.C_32_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_32_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_32_3")]; x <- cbind(x, "Intensity.C_32_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_4_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_4_2")]; x <- cbind(x, "Intensity.C_4_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_4_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_4_3")]; x <- cbind(x, "Intensity.C_4_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_4_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_4_4")]; x <- cbind(x, "Intensity.C_4_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_64_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_64_1")]; x <- cbind(x, "Intensity.C_64_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_64_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_64_3")]; x <- cbind(x, "Intensity.C_64_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_64_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_64_4")]; x <- cbind(x, "Intensity.C_64_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_64_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_64_5")]; x <- cbind(x, "Intensity.C_64_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]C_8_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]C_8_1")]; x <- cbind(x, "Intensity.C_8_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_0_5_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_0_5_1")]; x <- cbind(x, "Intensity.R_0_5_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_0_5_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_0_5_2")]; x <- cbind(x, "Intensity.R_0_5_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_0_5_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_0_5_5")]; x <- cbind(x, "Intensity.R_0_5_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_1_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_1_2")]; x <- cbind(x, "Intensity.R_1_2" = tmp)
x <- x[, !(colnames(x) %in% c("Intensity.R_1_3gap"))] ## remove
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_1_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_1_5")]; x <- cbind(x, "Intensity.R_1_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_16_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_16_1")]; x <- cbind(x, "Intensity.R_16_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_16_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_16_3")]; x <- cbind(x, "Intensity.R_16_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_2_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_2_5")]; x <- cbind(x, "Intensity.R_2_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_32_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_32_2")]; x <- cbind(x, "Intensity.R_32_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_32_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_32_3")]; x <- cbind(x, "Intensity.R_32_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_32_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_32_4")]; x <- cbind(x, "Intensity.R_32_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_32_5")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_32_5")]; x <- cbind(x, "Intensity.R_32_5" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_4_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_4_2")]; x <- cbind(x, "Intensity.R_4_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_4_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_4_3")]; x <- cbind(x, "Intensity.R_4_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_64_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_64_1")]; x <- cbind(x, "Intensity.R_64_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_64_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_64_3")]; x <- cbind(x, "Intensity.R_64_3" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_64_4")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_64_4")]; x <- cbind(x, "Intensity.R_64_4" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_8_1")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_8_1")]; x <- cbind(x, "Intensity.R_8_1" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_8_2")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_8_2")]; x <- cbind(x, "Intensity.R_8_2" = tmp)
tmp <- apply(x[, grep(colnames(x), pattern = "Intensity[.]R_8_3")], 1, mean, na.rm = TRUE)
x <- x[,-grep(colnames(x), pattern = "Intensity[.]R_8_3")]; x <- cbind(x, "Intensity.R_8_3" = tmp)
x <- x[, !(colnames(x) %in% c("Intensity"))] ## remove
x <- x[, !(colnames(x) %in% c("Intensity.0_6", "Intensity.0_7", "Intensity.0_8"))]
colnames(x) <- gsub(colnames(x), pattern = "C_0_5_", replacement = "C_05_") ## rename columns
colnames(x) <- gsub(colnames(x), pattern = "R_0_5_", replacement = "R_05_") ## rename columns

## remove first line since it contains the time
x <- x[-1, ]

## remove outlier Intensity.0_5
x <- x[, !(colnames(x) %in% c("Intensity.0_5"))]

x_old <- x

x <- apply(x, 2, function(x) x / quantile(x, .75, na.rm=T))

## take log values
#x <- log2(x+1)
## remove rows that contain only NaN values 


ind_batch1 <- grep(colnames(x), pattern = "0_1|05_1|1_1|2_1|4_1|8_1|16_1|32_1|64_1")
ind_batch2 <- grep(colnames(x), pattern = "0_2|05_2|1_2|2_2|4_2|8_2|16_2|32_2|64_2")
ind_batch3 <- grep(colnames(x), pattern = "0_3|05_3|1_3|2_3|4_3|8_3|16_3|32_3|64_3")
ind_batch4 <- grep(colnames(x), pattern = "0_4|05_4|1_4|2_4|4_4|8_4|16_4|32_4|64_4")
ind_batch5 <- grep(colnames(x), pattern = "0_5|05_5|1_5|2_5|4_5|8_5|16_5|32_5|64_5")

## remove the rows that have more than 4 missing values per row
ind_remove <- as.numeric(which(apply(x[, ind_batch1], 1, function(x) sum(is.na(x))) > 5))
ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch2], 1, function(x) sum(is.na(x))) > 5)))
ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch3], 1, function(x) sum(is.na(x))) > 5)))
ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch4], 1, function(x) sum(is.na(x))) > 5)))
ind_remove <- c(ind_remove, as.numeric(which(apply(x[, ind_batch5], 1, function(x) sum(is.na(x))) > 5)))
ind_remove <- unique(ind_remove)
x <- x[-ind_remove, ]

ind_batch_l <- list(ind_batch1, ind_batch2, ind_batch3, ind_batch4, ind_batch5)


################################################################################
## Day Normalization
x_DNtemp <- x

for (i in 1:length(ind_batch_l)) {
    subdata <- x[, ind_batch_l[[i]] ]
    metab.med <- apply(subdata,1,median,na.rm=T)
    x_DNtemp[, ind_batch_l[[i]] ] <- sweep(subdata,1,metab.med , FUN="-")
}

## add median per metabolite of whole experiment.
x.med.all <- apply(x ,1, median, na.rm=T)
x_DN <- x_DNtemp + x.med.all

################################################################################
## sample median normalization
sample.med <- apply(x_DN,2, median, na.rm=T)
x_DMN <- sweep(x_DN,2, sample.med, FUN="-") + median(sample.med)


## PCA
library(pcaMethods)
x <- as.matrix(x)
x <- log2(x+1)

x_DN <- as.matrix(x_DN)
x_DMN <- as.matrix(x_DMN)

stopifnot(is.numeric(x_DN))
stopifnot(is.numeric(x_DMN))

x <- ifelse(is.na(x), NA, x)
x_DN <- ifelse(is.na(x_DN), NA, x_DN)
x_DMN <- ifelse(is.na(x_DMN), NA, x_DMN)


x_qu <- sweep(x, MARGIN = 2, FUN = "/", STATS = apply(x, 2, function(y) quantile(y, 0.75, na.rm = T)) )

rpcaNN <- pca(t(x), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455, na.rm = TRUE)
rpcaDN <- pca(t(x_DN), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)
rpcaDMN <- pca(t(x_DMN), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)
rpcaQua <- pca(t(x_qu), method="ppca", center=TRUE, scale='none', nPcs=5, seed=3455)

# define colors (experiment dependent)
cls <- unlist(lapply(strsplit(colnames(x), split = "_"), function(x) paste(x[1], x[2], sep = "_")) )
cls[grep(cls, pattern = "Intensity[.]0_")] <- "Intensity_0"
cls <- unlist(lapply(strsplit(cls, split = "_"), "[", 2))
cls <- as.factor(cls)
cls <- as.numeric(cls)

#pchCR
pchCR <- numeric(length(cls))
pchCR[grep(colnames(x), pattern = "Intensity[.]C")] <- 1
pchCR[grep(colnames(x), pattern = "Intensity[.]R")] <- 2
pchCR[grep(colnames(x), pattern = "Intensity[.]0")] <- 3

#pdf('PCA_Mix_sample_NN.pdf',7,7)
slplot(rpcaNN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
slplot(rpcaNN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
slplot(rpcaNN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
#dev.off()

#pdf('PCA_Mix_sample_DN.pdf',7,7)
slplot(rpcaDN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
slplot(rpcaDN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
slplot(rpcaDN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2 )
#dev.off()


#pdf('PCA_Mix_sample_DMN.pdf',7,7)
slplot(rpcaDMN, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
slplot(rpcaDMN, pcs=c(1,3), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
slplot(rpcaDMN, pcs=c(2,3), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )
#dev.off()

slplot(rpcaQua, pcs=c(1,2), scoresLoadings = c(T, F), scol=cls, scex=1.2, sl = NULL, spch = pchCR )



