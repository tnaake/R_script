## load required libraries
library(igraph)
library(MSnbase)

################################### functions ##################################
#' @name shiftMatrix
#' @title Shift columns of a matrix by n and set added columns to def
#' @description \code{shiftMatrix} shifts columns of a matrix by \code{n} and 
#' sets the added columns to \code{def}.
#' @usage shiftMatrix(mat, x, n, def=NA)
#' @param mat matrix
#' @param x numeric, col indices to shift
#' @param n numeric
#' @param def character/numeric, replacement value for added columns
#' @details Helper function for \code{matchSpectra}
#' @return matrix with all combinations of shifted rows
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' shiftMatrix(mat, x, n, def=NA)
shiftMatrix <- function(mat, x, n, def=NA){
    if (n==0) { res <- mat }
    if (n<0) {
        n <- abs(n)
        res <- mat[,x[seq_len(length(x)-n)+n] ]##, rep(def, n))
        res <- cbind(matrix(res, nrow=nrow(mat), byrow=FALSE), 
                      matrix(def, nrow=nrow(mat), ncol=n))
    } else {
        res <- mat[,x[seq_len(length(x)-n)]]
        res <- cbind(matrix(def, nrow=nrow(mat), ncol=n), 
                      matrix(res, nrow(mat), byrow=FALSE))
    }
    return(res)
}

## test_that
library("testthat")
mat_l <- matrix(letters[1:18], ncol=6, nrow=3)
## tests
shiftMatrix(mat_l, x=c(2, 4, 6), n=-1)
matrix(c("j", "p", NA, "k", "q", NA, "l", "r", NA), ncol=3, nrow=3, byrow=TRUE)
shiftMatrix(mat_l, x=c(2, 4, 6), n=-2)
matrix(c("p", NA, NA, "q", NA, NA, "r", NA, NA), ncol=3, nrow=3, byrow=TRUE)

shiftMatrix(mat_l, x=c(2, 4, 6), n=1)
matrix(c(NA, "d", "j", NA, "e", "k", NA, "f", "l"), ncol=3, nrow=3, byrow=TRUE)
shiftMatrix(mat_l, x=c(2, 4, 6), n=2)
matrix(c(NA, NA, "d", NA, NA, "e", NA, NA, "f"), ncol=3, nrow=3, byrow=TRUE)



## taken from MetCirc
#' @name compare_Spectra
#' @title Create similarity matrix from MSnbase::Spectra object
#' @description compare_Spectra creates a similarity matrix of all 
#' Spectrum objects in \code{object}
#' @usage compare_Spectra(object, fun, ...)
#' @param object \code{Spectra}
#' @param fun function or character, see ?MSnbase::compareSpectra for further
#' information
#' @param ... further arguments passed to compareSpectra
#' @details Function inspired by compareSpectra.OnDiskMSnExp. Possibly 
#' transfer to MSnbase.
#' @return matrix with pair-wise similarities
#' @author Thomas Naake (inspired by compareSpectra.OnDiskMSnExp)
#' @examples 
#' data("spectra", package="MetCirc")
#' compare_Spectra(spectra_tissue[1:10], fun="dotproduct")
#' @export
compare_Spectra <- function(object, fun, ...) {
  
  nm <- names(object)
  cb <- combn(nm, 2)
  cb <- apply(cb, 2, function(x) compareSpectra(object[[x[1]]], object[[x[[2]]]], fun=fun, ...)) ## "dotproduct"
  
  m <- matrix(NA, length(object), length(object),
              dimnames=list(nm, nm))
  ## fill lower triangle of the matrix
  m[lower.tri(m)] <- cb
  ## copy to upper triangle
  for (i in 1:nrow(m)) {
    m[i, ] <- m[, i]
  }
  
  diag(m) <- 1
  
  return(m)
}


## taken from MetCirc
#' @name normalizeddotproduct
#' @title Calculate the normalized dot product
#' @description Calculate the normalized dot product (NDP)
#' @usage normalizeddotproduct(x, y, m=0.5, n=2, ...)
#' @param x \code{Spectrum2} object from \code{MSnbase} containing
#' intensity and m/z values, first MS/MS spectrum
#' @param y \code{Spectrum2} object from \code{MSnbase} containing 
#' intensity and m/z values, second MS/MS spectrum
#' @param m \code{numeric}, exponent to calculate peak intensity-based weights
#' @param n \code{numeric}, exponent to calculate m/z-based weights
#' @param ... further arguments passed to MSnbase:::bin_Spectra
#' @details The normalized dot product is calculated according to the 
#' following formula: 
#'  \deqn{NDP = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#'  with \eqn{W = [ peak intensity] ^{m} \cdot [m/z]^n}. For further information 
#'  see Li et al. (2015): Navigating natural variation in herbivory-induced
#'  secondary metabolism in coyote tobacco populations using MS/MS structural 
#'  analysis. PNAS, E4147--E4155. \code{normalizeddotproduct} returns a numeric 
#'  value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical. Arguments can be passed to 
#' the function \code{MSnbase:::bin_Spectra}, e.g. to set the width of bins
#' (binSize). 
#' Prior to calculating \deqn{W_{S1}} or \deqn{W_{S2}}, all intensity values 
#' are divided by the maximum intensity value. 
#' @return normalizeddotproduct returns a numeric similarity coefficient between 0 and 1
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' x <- spectra_tissue[[1]]
#' y <- spectra_tissue[[2]]
#' normalizeddotproduct(x, y, m=0.5, n=2, binSize=0.01) 
#' @export
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

#' @name matchSpectra
#' @title Match two spectra using bipartite networks and combinatorics
#' @description \code{matchSpectra} takes two objects, \code{x} and 
#' \code{y} as input that contain spectral information. The matching 
#' is a multi-step procedure: 
#' 1) filtering based on \code{ppm},
#' 2) retain order of matches between features of \code{x} and 
#' \code{y} (remove crossing edges that violate the order of matching
#' m/z),
#' 3) calculate all combinations of the remaining possibilities.  
#' @usage matchSpectra(x, y, ppm=20, fun=normalizeddotproduct, ...)
#' @param x matrix, the first row contains m/z value and the second 
#' row contains the corresponding intensity values                            --> x
#' @param  y matrix, the first row contains m/z value and the second 
#' row contains the corresponding intensity values                          --> y
#' @param ppm numeric, tolerance parameter in ppm to match corresponding peaks
#' @param fun function
#' @param ... additional parameters passed to compare_Spectra
#' @details Objective function is highest similarites between the two 
#' spectral objects, i.e. \code{fun} is calculated over all combinations and 
#' the similarity of the combination that yields the highest similarity is 
#' returned. 
#' @return list with elements x and y each being a matrix (columns "mz" 
#' each row (peak) in x matching the row (peak) in y
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' matchSpectra(x, y, ppm=20, fun=normalizeddotproduct, ...)
## function to match spectrum objects using bipartite networks and combinations
## to match peaks --> objective function is highest similarity between y
matchSpectra <- function(x, y, ppm=20, fun=normalizeddotproduct, ...) {
    
    if (!is.matrix(x)) stop("x is not a matrix")
    if (!is.matrix(y)) stop("y is not a matrix")
    if (!all(colnames(x) %in% c("mz", "intensity"))) stop("colnames(x) are not 'mz' and 'intensity'")
    if (!all(colnames(y) %in% c("mz", "intensity")))stop("colnames(y) are not 'mz' and 'intensity'")
    
    ## re-set colnames for x and y and order
    rownames(x) <- paste("sp1_", 1:nrow(x),sep="")
    rownames(y) <- paste("sp2_", 1:nrow(y),sep="")
    if (nrow(x) > 1) x <- x[order(x[, 1]), ]
    if (nrow(y) > 1) y <- y[order(y[, 1]), ]
    
    ## 1) create adjacency matrix and remove edges within x and 
    ## within y
    w <- matrix(1, ncol=nrow(x)+nrow(y), 
        nrow=nrow(x)+nrow(y),
        dimnames=list(c(rownames(x), rownames(y)),
            c(rownames(x), rownames(y))
    ))
    
    w[rownames(x), rownames(x)] <- 0
    w[rownames(y), rownames(y)] <- 0
    
    ## 2) remove edges that are not in a certain range
    ppm_1_1 <- x[,1] / abs(ppm / 10 ^ 6  - 1 ); names(ppm_1_1) <- rownames(x) 
    ppm_1_2 <- x[,1] / abs(ppm / 10 ^ 6  + 1 ); names(ppm_1_2) <- rownames(x)
    ppm_2_1 <- y[,1] / abs(ppm / 10 ^ 6  - 1 ); names(ppm_2_1) <- rownames(y) 
    ppm_2_2 <- y[,1] / abs(ppm / 10 ^ 6  + 1 ); names(ppm_2_2) <- rownames(y)

    mat1 <- apply(as.matrix(ppm_1_2), 1, function(a) a <= c(ppm_1_1, ppm_2_1))
    mat2 <- apply(as.matrix(ppm_1_1), 1, function(a) a >= c(ppm_1_2, ppm_2_2))

    link_ppm <- mat1*mat2
    w[rownames(link_ppm), colnames(link_ppm)] <- link_ppm
    w[colnames(link_ppm), rownames(link_ppm)] <- t(link_ppm)
    w[rownames(x), rownames(x)] <- 0
    w[rownames(y), rownames(y)] <- 0

    ## obtain network components from w
    net <- graph_from_adjacency_matrix(w, weighted=NULL, mode="undirected")
    comp <- components(net)
    
    ## 3) get all possible combinations within one component
    ## res will contain all combinations within one component
    res <- vector("list", length(comp$csize)) 
    
    ## write to res combinations where csize == 2
    inds_1 <- which(comp$csize == 1)
    res[inds_1] <- lapply(inds_1, function(a) {
        ms <- names(which(comp$membership == a))
        if (grepl(x=ms, pattern="sp1")) {ms <- c(ms, NA) } else{ms <- c(NA, ms)}
        matrix(ms, ncol=2)
    })
    ## write to res combinations where csize == 2
    inds_2 <- which(comp$csize == 2)
    res[inds_2] <- lapply(inds_2, function(a) matrix(names(which(comp$membership == a)), ncol=2))

    ## get combinations where csize > 2
    inds <- which(comp$csize > 2)
    ####################################################################################################
    for (i in inds) {
        
        ## separate component and create two matrices from x and 
        ## y that only contain component menbers
        names_ind_i <- names(which(comp$membership == i))
        x_ind <- x[names_ind_i[names_ind_i %in% rownames(x)], ]
        x_ind <- matrix(x_ind, ncol=2,
            dimnames=list(names_ind_i[names_ind_i %in% rownames(x)], c("", "")))
        y_ind <- y[names_ind_i[names_ind_i %in% rownames(y)], ]
        y_ind <- matrix(y_ind, ncol=2, 
            dimnames=list(names_ind_i[names_ind_i %in% rownames(y)], c("", ""))) 
        
        ## allocate to c1 and c2 the names (colnames of x and y 
        ## that are in the specific component)
        if (nrow(x_ind) < nrow(y_ind)) {
            c1 <- rownames(x_ind); c2 <- rownames(y_ind)
            c1 <- c(c1, rep("NA", length(c2)-length(c1))) 
            c1_c2 <- lapply(c1, function(a) expand.grid(a, c2))
        } else {
            c1 <- rownames(y_ind); c2 <- rownames(x_ind)
            c1 <- c(c1, rep("NA", length(c2)-length(c1))) 
            c1_c2 <- lapply(c2, function(a) expand.grid(a, c1))
        }
        ##c1 <- c(c1, rep("NA", length(c2)-length(c1))) #################################################################
        ##c1_c2 <- lapply(c1, function(a) expand.grid(a, c2))
        c1_c2_paste <- lapply(c1_c2, function(a) apply(a, 1, function(b) paste(b, collapse=" & ")))
    
        ## calculate all possible combinations
        res_i <- as.matrix(expand.grid(c1_c2_paste, stringsAsFactors=FALSE))
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
        if (ncol(res_i) > 2) { ## do this if there is only ############################################
            ## a mapping multiple to multiple
        
            ## check order of sp2s, they have to ascend, remove those that
            ## do not ascend
            crosses <- lapply(as.data.frame(t(res_i[, seqs]), stringsAsFactors=FALSE), 
                function(a) as.numeric(substr(a, 5, nchar(a))))
            ## check where all are TRUE, remove before NAs
            crosses <- lapply(crosses, function(a) {
              a <- a[!is.na(a)]
              all(a == sort(a))
            })
            ##crosses <- lapply(seq_along(crosses[-1]), function(a) crosses[[a]] <= crosses[[a+1]])
            res_i <- matrix(res_i[unlist(crosses), ], ncol=ncol(res_i))
            
            ## create shifted matrices, do not use last element since it contains
            ## only NAs 
            shift_right <- lapply(seqs[-length(seqs)], function(a) shiftMatrix(res_i, seqs, a/2))
            shift_left <- lapply(seqs[-length(seqs)], function(a) shiftMatrix(res_i, seqs, -a/2))
    
            ## write to matrix
            mat_shift_left <- lapply(shift_left, function(a) {
                res_i[, seqs] <- a
                return(res_i)
            })
            mat_shift_right <- lapply(shift_right, function(a) {
                res_i[, seqs] <- a
                return(res_i)
            })
            
            ## rbind the lists
            mat_shift_left <- do.call(rbind, mat_shift_left)
            mat_shift_right <- do.call(rbind, mat_shift_right)
      
            ## the ones that are shifted out, link to NA and bind to NA
            ## for mat_shift_left (take from shift_right and inverse)
            add_left <- shift_right[length(shift_right):1]
            ##add_left <- lapply(add_left, function(x) x[nrow(x):1],)
            ##add_left <- lapply(add_left, function(x) x[length(x):1])
            add_left <- lapply(add_left, function(x) x[,-1])
            
            mat_add_left <- mat_add_right <- matrix("NA", ncol=length(shift_left)*2, nrow=nrow(mat_shift_left))
            mat_add_left[, seqs[-length(seqs)]] <- unlist(add_left)
            
            ## for mat_shift_right (take from shift_left and inverse)
            add_right <- shift_left[length(shift_left):1]
            add_right <- lapply(add_right, function(x) matrix(x[,ncol(x):1], ncol=ncol(x)))
            add_right <- lapply(add_right, function(x) x[,-1])
            mat_add_right[, seqs[-length(seqs)]] <- unlist(add_right)
            
            ## cbind with mat_add_left and mat_add_right
            mat_shift_left <- cbind(mat_shift_left, mat_add_left)
            mat_shift_right <- cbind(mat_shift_right, mat_add_right)
            
            #apply(mat_shift_left, 1, table)
            #apply(mat_shift_right, 1, table)
            
            ## add to res_i
            res_i <- cbind(res_i, matrix("NA", nrow=nrow(res_i), ncol=ncol(mat_add_left)))
            
            ## assign to res
            res[[i]] <- rbind(res_i, mat_shift_left, mat_shift_right)
        } else {
            res[[i]] <- res_i
        }
    }

    ## create combinations between rows 
    res_paste <- lapply(res, function(a) apply(a, 1, paste, collapse=" & "))

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
    sim <- apply(res_exp, 1, function(a) {
        sp1_ind <- a[seq(1, ncol(res_exp), 2)] 
        sp2_ind <- a[seq(2, ncol(res_exp), 2)] 
        
        ## remove elements when two "NA" are at the same position 
        ## (they were created artificially in the above step when pushing 
        ## sp1/sp2 out and cbinding the pushed out from the matrix)
        ind_remove <- sp1_ind == "NA" & sp2_ind == "NA"
        sp1_ind <- sp1_ind[!ind_remove]
        sp2_ind <- sp2_ind[!ind_remove]
        
        ## create vectors that store mz and intensity of combination
        mz1 <- rep(NA, length(sp1_ind))
        int1 <- numeric(length(sp1_ind))
        mz2 <- rep(NA, length(sp2_ind))
        int2 <- numeric(length(sp2_ind))
        mz1[sp1_ind != "NA"] <- x[sp1_ind[sp1_ind != "NA"], "mz"]
        mz2[sp2_ind != "NA"] <- y[sp2_ind[sp2_ind != "NA"], "mz"]
        int1[sp1_ind != "NA"] <- x[sp1_ind[sp1_ind != "NA"], "intensity"]
        int2[sp2_ind != "NA"] <- y[sp2_ind[sp2_ind != "NA"], "intensity"]
        
        mz1_old <- mz1
        mz2_old <- mz2
        ## replace NA value in mz1 and mz2 by average of neighbours in order to 
        ## secure that mz are strictly monotonically increasing
        ## PROBLEM: when any other function is used this can result to errors
        ## PROBLEM: what is better approach? take mz from the other spectrum at that
        ## position --> could be smaller or greater than neighbour m/z values
        ## --> when creating y mz values are ordered increasingly
  
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
    
        ## problem y reorders mz and associated ntensities!!!!
        
        ## create y objects
        S2_1 <- new("Spectrum2", precursorMz=max(mz1, na.rm=TRUE), mz=mz1, intensity=int1)
        S2_2 <- new("Spectrum2", precursorMz=max(mz2, na.rm=TRUE), mz=mz2, intensity=int2)
        ## calculate similarity
        value <- compareSpectra(S2_1, S2_2, fun=fun, binning=FALSE, ...)##compareSpectra(S2_1, S2_2, fun=fun, binning=FALSE, ...)
        
        ## cbind spectra
        x <- cbind(mz=mz1_old, intensity=int1)
        y <- cbind(mz=mz2_old, intensity=int2)
        
        l <- list(value=value, x=x, y=y)
        return(l)
    })
    sim_value <- unlist(lapply(sim, "[[", "value"))
    sim_max <- sim[[which.max(sim_value)]]
    
    ## sort according to ascending order
    x_max <- sim_max[["x"]]
    y_max <- sim_max[["y"]]
    
    sort_x_y <- lapply(seq_len(nrow(x_max)), function(a) 
      paste(sort(c(x_max[a, "mz"], y_max[a, "mz"])), collapse="_")
    )
    sort_x_y <- order(unlist(sort_x_y))
    x_max <- x_max[sort_x_y,]
    y_max <- y_max[sort_x_y,]
    
    l <- list(x=x_max, y=y_max)
    return(l)

}

################################### workflow ###################################
## two example spectras
spectrum1 <- matrix(c(c(100, 200, 200.001, 400, 400.00005, 400.0001, 400.00011, 400.00012),
                      c(1, 1, 1, 1, 1.5, 1.2, 1.0, 1.0)), ncol=2, nrow=8, byrow=FALSE)
colnames(spectrum1) <- c("mz", "intensity")

spectrum2 <- matrix(c(c(100.001, 199.999, 200.0005, 399.99998, 399.999999, 400.00005, 400.00006),
                      c(1, 1, 1, 1.5, 2, 2.5, 2.4)), ncol=2, nrow=7, byrow=FALSE)
colnames(spectrum2) <- c("mz", "intensity")

## should give
##100.001 <-> 100
##100.002 <-> NA --> gives higher score
##NA      <-> 200.0
##300.01  <-> 300.002
##300.02  <-> 300.0255 --> gives higher score
##NA      <-> 300.0250
matchSpectra(x=spectrum1[-(2:3),], y=spectrum2, n=1, m=0.2) ## 0.998162

sp1 <- t(as.matrix(spectrum1[-c(2:4),]))
sp2 <- spectrum2[-1, ]



## unit tests via test_that
library("testthat")
## create example spectrum1 and spectrum2 and perform tests
spectrum1 <- matrix(c(c(100.001, 100.002, 300.01, 300.02),
                      c(1, 1, 1, 1)), ncol=2, nrow=4, byrow=FALSE)
colnames(spectrum1) <- c("mz", "intensity")

spectrum2 <- matrix(c(c(100.0, 200.0, 300.002, 300.025, 300.0255),
                      c(1, 1, 1, 1, 1)), ncol=2, nrow=5, byrow=FALSE)
colnames(spectrum2) <- c("mz", "intensity")

spectrum1_match <- matrix(c(100.002, 100.001, NA, 300.01, 300.02, 
    NA, 1, 1, 0, 1, 1, 0), ncol=2, nrow=6, byrow=FALSE, dimnames=list(NULL, c("mz", "intensity")))
spectrum2_match <- matrix(c(100.0, NA, 200.0, 300.002, 300.025, 300.0255,
    1, 0, 1, 1, 1, 1), ncol=2, nrow=6, byrow=FALSE, dimnames=list(NULL, c("mz", "intensity")))

test_that("", {
  expect_equal(matchSpectra(x=spectrum1, y=spectrum1), list(x=spectrum1, y=spectrum1))
  expect_equal(matchSpectra(x=spectrum2, y=spectrum2), list(x=spectrum2, y=spectrum2))
  expect_equal(matchSpectra(x=spectrum1, y=spectrum2), list(x=spectrum1_match,
                                                            y=spectrum2_match))
  expect_error(matchSpectra(x=spectrum1[1,], y=spectrum2))
  expect_error(matchSpectra(x=spectrum1))
  expect_error(matchSpectra(y=spectrum2))
  expect_error(matchSpectra(x=spectrum1, y=spectrum2, fun=max))
})
