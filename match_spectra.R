## load required libraries
library(igraph)
library(MSnbase)

## two example spectras
spectrum1 <- matrix(c(c(100, 200, 200.001, 400, 400.00005, 400.0001),
                      c(1, 1, 1, 1, 1.5, 1.2)), ncol=6, nrow=2, byrow=TRUE)

spectrum2 <- matrix(c(c(100.001, 199.999, 200.0005, 399.99998, 399.999999, 400.00005, 400.00006),
                      c(1, 1, 1, 1.5, 2, 2.5, 2.4)), ncol=7, nrow=2, byrow=TRUE)



## helper function: shift columns of matrix by n and set to def
shift_matrix <- function(mat, x, n, def=NA){
    if(n==0){ res <- mat}
    if(n<0){
        n <- abs(n)
        res <- mat[,x[seq_len(length(x)-n)+n] ]##, rep(def, n))
        res <- cbind(res, matrix(NA, nrow=nrow(mat), ncol=n))
    }else{
        ##res <- c(rep(def, n), 
        res <- mat[,x[seq_len(length(x)-n)]]
        res <- cbind(matrix(NA, nrow=nrow(mat), ncol=n), res)
    }
    return(res)
}

## function to calculate dotproduct
normalizeddotproduct <- function(x, y, m=0.5, n=2, binning=FALSE, ...) {
    ## normalize to % intensity
    x@intensity <- x@intensity / max(x@intensity)*100
    y@intensity <- y@intensity / max(y@intensity)*100
    
    if (binning) {
        binnedSpectra <- MSnbase:::bin_Spectra(x, y, ...)   
        inten <- lapply(binnedSpectra, intensity)
        mz <- lapply(binnedSpectra, mz)    
    } else {
        inten <- list(x@intensity, y@intensity)
        mz <- list(x@mz, y@mz)
    }
    
    
    ws1 <- inten[[1]] ^ m * mz[[1]] ^ n
    ws2 <- inten[[2]] ^ m * mz[[2]] ^ n
    
    sum( ws1*ws2, na.rm=TRUE) ^ 2 / ( sum( ws1^2, na.rm=TRUE) * sum( ws2^2, na.rm=TRUE ) )
}

## function to match Spectrum2 objects using bipartite networks and combinations
## to match peaks --> objective function is highest similarity between Spectrum2
matchSpectrum2 <- function(spectrum1, spectrum2, ppm=20, fun=normalizeddotproduct, ...) {
    
    if (!is.matrix(spectrum1)) stop("spectrum1 is not a matrix")
    if (!is.matrix(spectrum2)) stop("spectrum2 is not a matrix")
    
    ## re-set colnames for spectrum1 and spectrum2 and order
    colnames(spectrum1) <- paste("sp1_", 1:ncol(spectrum1),sep="")
    colnames(spectrum2) <- paste("sp2_", 1:ncol(spectrum2),sep="")
    if (ncol(spectrum1) > 1) spectrum1 <- spectrum1[, order(spectrum1[1, ])]
    if (ncol(spectrum2) > 1) spectrum2 <- spectrum2[, order(spectrum2[1, ])]
    
    ## 1) create adjacency matrix and remove edges within spectrum1 and 
    ## within spectrum2
    w <- matrix(1, ncol=ncol(spectrum1)+ncol(spectrum2), 
        nrow=ncol(spectrum1)+ncol(spectrum2),
        dimnames=list(c(colnames(spectrum1), colnames(spectrum2)),
            c(colnames(spectrum1), colnames(spectrum2))
    ))
    
    w[colnames(spectrum1), colnames(spectrum1)] <- 0
    w[colnames(spectrum2), colnames(spectrum2)] <- 0
    
    ## 2) remove edges that are not in a certain range
    ppm_1_1 <- spectrum1[1,] / abs(ppm / 10 ^ 6  - 1 ) 
    ppm_1_2 <- spectrum1[1,] / abs(ppm / 10 ^ 6  + 1 ) 
    ppm_2_1 <- spectrum2[1,] / abs(ppm / 10 ^ 6  - 1 ) 
    ppm_2_2 <- spectrum2[1,] / abs(ppm / 10 ^ 6  + 1 ) 

    mat1 <- apply(as.matrix(ppm_1_2), 1, function(x) x <= c(ppm_1_1, ppm_2_1))
    mat2 <- apply(as.matrix(ppm_1_1), 1, function(x) x >= c(ppm_1_2, ppm_2_2))

    link_ppm <- mat1*mat2
    w[rownames(link_ppm), colnames(link_ppm)] <- link_ppm
    w[colnames(link_ppm), rownames(link_ppm)] <- t(link_ppm)
    w[colnames(spectrum1), colnames(spectrum1)] <- 0
    w[colnames(spectrum2), colnames(spectrum2)] <- 0

    ## match order and remove 
    comp <- components(graph_from_adjacency_matrix(w))
    
    ## 3) get all possible combinations within one component
    ## res will contain all combinations within one component
    res <- vector("list", length(comp$csize)) 
    
    ## write to res combinations where csize == 2
    inds_1 <- which(comp$csize == 1)
    res[inds_1] <- lapply(inds_1, function(x) {
        ms <- names(which(comp$membership == x))
        if (grepl(x=ms, pattern="sp1")) {ms <- c(ms, NA) } else{ms <- c(NA, ms)}
        matrix(ms, ncol=2)
    })
    ## write to res combinations where csize == 2
    inds_2 <- which(comp$csize == 2)
    res[inds_2] <- lapply(inds_2, function(x) matrix(names(which(comp$membership == x)), ncol=2))

    ## get combinations where csize > 2
    inds <- which(comp$csize > 2)
    for (i in inds) {
        
        ## separate component and create two matrices from spectrum1 and 
        ## spectrum2 that only contain component menbers
        names_ind_i <- names(which(comp$membership == i))
        s1_ind <- spectrum1[, names_ind_i[names_ind_i %in% colnames(spectrum1)] ]
        s1_ind <- matrix(s1_ind, nrow=2,
            dimnames=list(c("", ""), names_ind_i[names_ind_i %in% colnames(spectrum1)]))
        s2_ind <- spectrum2[, names_ind_i[names_ind_i %in% colnames(spectrum2)] ]
        s2_ind <- matrix(s2_ind, nrow=2, 
            dimnames=list(c("", ""), names_ind_i[names_ind_i %in% colnames(spectrum2)]))
        

        ## allocate to x and y the names (colnames of spectrum1 and spectrum2 
        ## that are in the specific component)
        if (ncol(s1_ind) < ncol(s2_ind)) {
            x <- colnames(s1_ind); y <- colnames(s2_ind)
        } else {
            x <- colnames(s2_ind); y <- colnames(s1_ind)
        }
        x_y <- lapply(x, function(x) expand.grid(x, y))
        x_y_paste <- lapply(x_y, function(x) apply(x, 1, function(y) paste(sort(y), collapse=" & ")))
    
        ## calculate all possible combinations
        res_i <- as.matrix(expand.grid(x_y_paste, stringsAsFactors=FALSE))
        ## write rows to list
        res_i <- split(res_i, row(res_i))
        ## strsplit " & ", unlist and write to matrix 
        res_i <- lapply(res_i, strsplit, split=" & ")
        res_i <- lapply(res_i, unlist)
        res_i <- matrix(unlist(res_i), nrow=length(res_i), byrow=TRUE)
    
        ## remove rows that contain duplicated values 
        res_i <- res_i[!apply(apply(res_i, 1, duplicated), 2, any),]
        
        ## filtering for crossing matching: retain order of m/z
        seqs <- seq(2, ncol(res_i), by=2)
        if (nrow(res_i) > 2 & ncol(res_i) > 2) { ## do this if there is only 
            ## a mapping multiple to multiple
        
            crosses <- lapply(as.data.frame(res_i[, seqs], stringsAsFactors = F), 
                function(x) as.numeric(substr(x, 5, nchar(x))))
            crosses <- lapply(seq_along(crosses[-1]), function(x) crosses[[x]] <= crosses[[x+1]])
            ## check where all are TRUE: make a logical matrix and use apply(m, 2, all) 
            ## shift columns by +2, +4, +6, etc. and fill with NA
            crosses <-  matrix(unlist(crosses), nrow=length(crosses), byrow=TRUE)
            res_i <- res_i[apply(crosses, 2, all),]
            res_i <- matrix(res_i, ncol=length(seqs)*2)
    
            ## create shifted matrices, do not use last element since it contains
            ## only NAs 
            shift_right <- lapply(seqs[-length(seqs)], function(x) shift_matrix(res_i, seqs, x/2))
            shift_left <- lapply(seqs[-length(seqs)], function(x) shift_matrix(res_i, seqs, -x/2))
    
            ## write to matrix
            mat_shift_left <- lapply(shift_left, function(x) {
                res_i[, seqs] <- x
                return(res_i)
            })
            mat_shift_right <- lapply(shift_right, function(x) {
                res_i[, seqs] <- x
                return(res_i)
            })
    
            ## rbind the lists
            mat_shift_left <- do.call(rbind, mat_shift_left)
            mat_shift_right <- do.call(rbind, mat_shift_right)
    
            ## assign to res
            res[[i]] <- rbind(res_i, mat_shift_left, mat_shift_right)
        } else {
            res[[i]] <- res_i
        }
    }
    

    ## create combinations between rows 
    res_paste <- lapply(res, function(x) apply(x, 1, paste, collapse=" & "))

    ## 4) calculate all possible combinations between the components
    res_exp <- as.matrix(expand.grid(res_paste, stringsAsFactors=FALSE))
    ## write rows to list
    res_exp <- split(res_exp, row(res_exp))
    ## strsplit " & ", unlist and write to matrix 
    res_exp <- lapply(res_exp, strsplit, split=" & ")
    res_exp <- lapply(res_exp, unlist)
    res_exp <- matrix(unlist(res_exp), nrow=length(res_exp), byrow=TRUE)

    ## 5) go through every row and calculate score: row with highest score is 
    ## the best match
    sim <- apply(res_exp, 1, function(x) {
        sp1_ind <- x[seq(1, ncol(res_exp), 2)] 
        sp2_ind <- x[seq(2, ncol(res_exp), 2)] 
    
        mz1 <- rep(NA, length(sp1_ind))
        int1 <- numeric(length(sp1_ind))
        mz2 <- rep(NA, length(sp2_ind))
        int2 <- numeric(length(sp2_ind))
        mz1[sp1_ind != "NA"] <- spectrum1[1, sp1_ind[sp1_ind != "NA"]]
        mz2[sp2_ind != "NA"] <- spectrum2[1, sp2_ind[sp2_ind != "NA"]]
        int1[sp1_ind != "NA"] <- spectrum1[2, sp1_ind[sp1_ind != "NA"]]
        int2[sp2_ind != "NA"] <- spectrum2[2, sp2_ind[sp2_ind != "NA"]]

        ## replace NA value in mz1 and mz2 by average of neighbours in order to 
        ## secure that mz are strictly monotonically increasing
        ## PROBLEM: when any other function is used this can result to errors
        ## PROBLEM: what is better approach? take mz from the other spectrum at that
        ## position --> could be smaller or greater than neighbour m/z values
        ## --> when creating Spectrum2 mz values are ordered increasingly
  
        ## if the first element is NA set it to minimum value - 0.01 or 0
        if (is.na(mz1[1])) mz1[1] <- max(min(mz1, na.rm=TRUE) - 0.01, 0)
        if (is.na(mz2[1])) mz2[1] <- max(min(mz2, na.rm=TRUE) - 0.01, 0)
    
        ## if there are more NA values "impute" them by taking preceding 
        ## value and adding small value in order to keep the order
        ind_na_1 <- which(is.na(mz1))
        ind_na_2 <- which(is.na(mz2))
        if (length(ind_na_1) > 0) {
            for (i in ind_na_1) mz1[i] <- mz1[i-1]+1e-30   
        }
        if (length(ind_na_2) > 0) {
            for (i in ind_na_2) mz2[i] <- mz2[i-1]+1e-30    
        }
    
        ## problem Spectrum2 reorders mz and associated ntensities!!!!
        
        ## create Spectrum2 objects
        S2_1 <- new("Spectrum2", precursorMz=max(mz1, na.rm=TRUE), mz=mz1, intensity=int1)
        S2_2 <- new("Spectrum2", precursorMz=max(mz2, na.rm=TRUE), mz=mz2, intensity=int2)
        ## calculate similarity
        compareSpectra(S2_1, S2_2, fun=fun, binning=FALSE, ...)
    })

    return(max(sim))

}


matchSpectrum2(spectrum1=as.matrix(spectrum1[,-c(2:6)]), spectrum2=spectrum2, n=1, m=0.2)

sp1 <- as.matrix(spectrum1[, -c(2:6)])
sp2 <- spectrum2



## unit tests via test_that
## expect_equal
## expect_equal
library("testthat")

matchSpectrum2(spectrum1=spectrum1[,1:2], spectrum2=spectrum2[,-c(1:3)])

test_that("", {
  expect_equal(matchSpectrum2(spectrum1=spectrum1, spectrum2=spectrum1), 1)
  expect_equal(matchSpectrum2(spectrum1=spectrum1, spectrum2=spectrum2), 0.9957219)
  expect_equal(matchSpectrum2(spectrum1=spectrum1[,1:2], spectrum2=spectrum2[,-c(1:3)]), 0)
})