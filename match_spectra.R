## load required libraries
library(igraph)

################################### functions ##################################
#' @name shiftMatrix
#' @title Shift columns of a matrix by n and set added columns to def
#' @description `shiftMatrix` shifts columns of a matrix by `n` and 
#' sets the added columns to `def`.
#' @usage shiftMatrix(mat, x, n, def=NA)
#' @param mat `matrix`
#' @param x numeric, col indices to shift
#' @param n `numeric(1)`, gives the number by how many positions the columns 
#' `x` of `mat` are shifted
#' @param def `character(1)`/`numeric(1)`, replacement value for added columns
#' @details helper function for `graphPeaks`
#' @return matrix with all combinations of shifted rows, only returns the 
#' columns `x` of `mat`
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' mat <- matrix(letters[1:18], ncol=6, nrow=3)
#' x <- c(2, 4, 6)
#' shiftMatrix(mat=mat, x=x, n=1, def=NA)
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

## tests with test_that
library("testthat")

## create test matrices
mat_l <- matrix(letters[1:18], ncol=6, nrow=3)
## n: negative, p: positive
mat_n1 <- matrix(c("j", "p", NA, "k", "q", NA, "l", "r", NA), ncol=3, nrow=3, byrow=TRUE)
mat_n2 <- matrix(c("p", NA, NA, "q", NA, NA, "r", NA, NA), ncol=3, nrow=3, byrow=TRUE)
mat_p1 <- matrix(c(NA, "d", "j", NA, "e", "k", NA, "f", "l"), ncol=3, nrow=3, byrow=TRUE)
mat_p2 <- matrix(c(NA, NA, "d", NA, NA, "e", NA, NA, "f"), ncol=3, nrow=3, byrow=TRUE)


test_that("", {
    expect_equal(shiftMatrix(mat_l, x=c(2, 4, 6), n=-1), mat_n1)
    expect_equal(shiftMatrix(mat_l, x=c(2, 4, 6), n=-2), mat_n2)
    expect_equal(shiftMatrix(mat_l, x=c(2, 4, 6), n=1), mat_p1)
    expect_equal(shiftMatrix(mat_l, x=c(2, 4, 6), n=2), mat_p2)
    expect_error(shiftMatrix(x=c(2,4,6), n=1, def=NA))
    expect_error(shiftMatrix(mat=mat_l, n=1, def=NA))
    expect_error(shiftMatrix(mat=mat_l, x=c(2,4,6), def=NA))
    expect_error(shiftMatrix(mat=mat_l, x=c(2,4,6,8,10), n=1, def=NA))
})



#' @name normalizeddotproduct
#' @title Calculate the normalized dot product
#' @description Calculate the normalized dot product (NDP)
#' @usage normalizeddotproduct(x, y, m=0.5, n=2, ...)
#' @param x `list`/`data.frame` of length 2 with m/z (`"mz"`) and corresponding intensity 
#' values (`"intensity"`)
#' @param y `list`/`data.frame` of length 2 with m/z (`"mz"`) and corresponding intensity 
#' values (`"intensity"`)
#' @param m `numeric(1)`, exponent to calculate peak intensity-based weights
#' @param n `numeric(1)`, exponent to calculate m/z-based weights
#' @details The normalized dot product is calculated according to the 
#' following formula: 
#'  \deqn{NDP = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#'  with \eqn{W = [ peak intensity] ^{m} \cdot [m/z]^n}. For further information 
#'  see Li et al. (2015): Navigating natural variation in herbivory-induced
#'  secondary metabolism in coyote tobacco populations using MS/MS structural 
#'  analysis. PNAS, E4147--E4155. `normalizeddotproduct` returns a numeric 
#'  value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical. 
#' Prior to calculating \deqn{W_{S1}} or \deqn{W_{S2}}, all intensity values 
#' are divided by the maximum intensity value. 
#' @return `numeric(1)`, `normalizeddotproduct` returns a numeric similarity 
#' coefficient between 0 and 1
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' x <- spectra_tissue[[1]]
#' y <- spectra_tissue[[2]]
#' normalizeddotproduct(x, y, m=0.5, n=2, binSize=0.01) 
#' @export
normalizeddotproduct <- function(x, y, m=0.5, n=2) {
    
    ## retrieve m/z and intensity from x and y
    mz1 <- x$mz
    mz2 <- y$mz
    inten1 <- x$intensity
    inten2 <- y$intensity
    
    if (length(mz1) != length(mz2)) stop("length(mz1) not equal to length(mz2)")
    if (length(inten1) != length(mz2)) stop("length(mz1) not equal to length(mz2)")
    if (length(mz1) != length(inten1)) stop("length(mz1) not equal to length(inten1)")
    
    ## normalize to % intensity
    inten1 <- inten1 / max(inten1, na.rm=TRUE)*100
    inten2 <- inten2 / max(inten2, na.rm=TRUE)*100
    
    ws1 <- inten1 ^ m * mz1 ^ n
    ws2 <- inten2 ^ m * mz2 ^ n
    
    sum( ws1*ws2, na.rm=TRUE) ^ 2 / ( sum( ws1^2, na.rm=TRUE) * sum( ws2^2, na.rm=TRUE ) )
}

#' @name graphPeaks
#' @title Match two spectra using bipartite networks and combinatorics
#' @description `graphPeaks` takes two objects, `x` and 
#' `y` as input that contain spectral information. The matching 
#' is a multi-step procedure: 
#' 1) filtering based on `ppm`,
#' 2) retain order of matches between features of `x` and 
#' `y` (remove crossing edges that violate the order of matching
#' m/z),
#' 3) calculate all combinations of the remaining possibilities.  
#' @usage graphPeaks(x, y, ppm=20, fun=normalizeddotproduct, ...)
#' @param x `matrix`, the first column (`"mz"`) contains m/z value and the 
#' second column (`"intensity"`) contains the corresponding intensity values                            
#' @param y `matrix`, the first column (`"mz"`) contains m/z value and the 
#' second column (`"intensity"`) contains the corresponding intensity values                          
#' @param ppm numeric, tolerance parameter in ppm to match corresponding peaks
#' @param fun function to calculate similarity between spectra
#' @param ... additional parameters passed to `fun`
#' @details Objective function is highest similarites between the two 
#' spectral objects, i.e. `fun` is calculated over all combinations and 
#' the similarity of the combination that yields the highest similarity is 
#' returned. 
#' @return list with elements `x` and `y` each being a matrix with columns `"mz"` 
#' and `"intensity"`. Each row (peak) in `x` matches the row (peak) in `y`
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' graphPeaks(x, y, ppm=20, fun=normalizeddotproduct, ...)
## function to match spectrum objects using bipartite networks and combinations
## to match peaks --> objective function is highest similarity between y
graphPeaks <- function(x, y, ppm=20, fun=normalizeddotproduct, ...) {
    
    if (!is.matrix(x)) stop("x is not a matrix")
    if (!is.matrix(y)) stop("y is not a matrix")
    if (mode(x) != "numeric") stop("mode(x) is not 'numeric'")
    if (mode(y) != "numeric") stop("mode(y) is not 'numeric'")
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
        if (ncol(res_i) > 2) { 
            ## do, if there is only a multiple mapping
        
            ## check order of sp2s, they have to ascend, remove those that
            ## do not ascend
            crosses <- lapply(as.data.frame(t(res_i[, seqs]), stringsAsFactors=FALSE), 
                function(a) as.numeric(substr(a, 5, nchar(a))))
            ## check where all are TRUE, remove before NAs
            crosses <- lapply(crosses, function(a) {
              a <- a[!is.na(a)]
              all(a == sort(a))
            })
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
        
        sp1 <- data.frame(mz=mz1, intensity=int1)
        sp2 <- data.frame(mz=mz2, intensity=int2)
        
        value <- fun(sp1, sp2, ...)
        
        l <- list(value=value, x=sp1, y=sp2)
        return(l)
        
    })
    sim_value <- unlist(lapply(sim, "[[", "value"))
    sim_max <- sim[[which.max(sim_value)]]
    
    # ## sort according to ascending order
    x_max <- as.matrix(sim_max[["x"]])
    y_max <- as.matrix(sim_max[["y"]])
    
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
graphPeaks(x=spectrum1, y=spectrum2, fun=normalizeddotproduct, n=1, m=0) 

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
  expect_equal(graphPeaks(x=spectrum1, y=spectrum1), list(x=spectrum1, y=spectrum1))
  expect_equal(graphPeaks(x=spectrum2, y=spectrum2), list(x=spectrum2, y=spectrum2))
  expect_error(graphPeaks(x=spectrum1[1,], y=spectrum2))
  expect_error(graphPeaks(x=spectrum1))
  expect_error(graphPeaks(y=spectrum2))
  expect_error(graphPeaks(x=spectrum1, y=spectrum2, fun=max))
}


## chcek integration with functions
x <- new("Spectrum2", )
y <- new("Spectrum2")
.compare_spectra <- function(x, y = NULL, MAPFUN = joinPeaks, tolerance = 0,
                             ppm = 20, FUN = cor, ...) {
  x_idx <- seq_along(x)
  y_idx <- seq_along(y)
  
  nx <- length(x_idx)
  ny <- length(y_idx)
  
  mat <- matrix(NA_real_, nrow = nx, ncol = ny,
                dimnames = list(spectraNames(x), spectraNames(y)))
  
  ## Might need some tuning - bplapply?
  ## This code duplication may be overengineering.
  if (nx >= ny) {
    for (i in x_idx) {
      px <- peaks(x[i])[[1L]]
      for (j in y_idx) {
        peak_map <- MAPFUN(px, peaks(y[j])[[1L]],
                           tolerance = tolerance, ppm = ppm, ...)
        mat[i, j] <- FUN(peak_map[[1L]][, 2L], peak_map[[2L]][, 2L],
                         ...)
      }
    }
  } else {
    for (j in y_idx) {
      py <- peaks(y[j])[[1L]]
      for (i in x_idx) {
        peak_map <- MAPFUN(peaks(x[i])[[1]], py,
                           tolerance = tolerance, ppm = ppm, ...)
        mat[i, j] <- FUN(peak_map[[1L]][, 2L], peak_map[[2L]][, 2L],
                         ...)
      }
    }
  }
  mat
}

