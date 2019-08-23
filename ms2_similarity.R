################################## functions ###################################

## functions to create MS/MS features from corresponding files that were 
## collected under different collision energies

#' @name assemblySpectra
#' @title Deconvolute final spectra from different collision energies
#' @description assemblySpectra takes files from different collision energies
#' as input and assembles a final deconvoluted spectra from it, taking the 
#' highest intensity values from corresponding fragments across the different 
#' collision energies. 
#' @usage assemlySpectra(spectra, aln, cols_aln, sample_name, tol)
#' @param spectra list of matrices, each object is from MS-DIAL and contains
#' information on the MS/MS spectrum
#' @param aln matrix, alignment file from MS-DIAL
#' @param cols_aln character, colnames in aln that correspond to the objects in
#' spectra
#' @param sample_name character, add the name to the names of the returned list
#' (Alignment.ID, retention time and m/z of the precursor, taken from aln, will
#' be automatically added to the name)
#' @param tol numeric, tolerance parameter
#' @details assemblySpectra will first bin within each spectra (su)
#' @return list of assembled spectra. Each spectra is a matrix that in the 
#' first column contains m/z values and in the second column contains the 
#' corresponding intensities
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' sp <- list(MK_hcd30, MK_hcd40, MK_hcd50)
#' cols_aln <- c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50")
#' assemblySpectra(spectra=sp, aln=aln, cols_aln=cols_aln, sample_name="MK_HCD")
assemblySpectra <- function(spectra=list(MK_hcd30, MK_hcd40, MK_hcd50), aln=aln,
    cols_aln=c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD", tol=0.01) {
    
    
    ## create list of vector nrow(aln) that will store the final deconvoluted 
    ## spectra, names refer to the aligned spectra in aln
    res <- vector("list", nrow(aln))
    
    
    ## iterate through aln and assemble final deconvoluted spectra
    res <- lapply(seq_len(nrow(aln)), function(i) {
        
        ## retrieve PeakID from alignment object aln per row, get IDs from all
        ## collision energies and assign as a list to inds_ev
        inds_ev <- lapply(cols_aln, function(x) aln[i, x])
        
        ## retrieve the corresponding row in spectra (for each collision energy)
        ## the entry will store information on the spectrum, use the "PeakID"
        ## as the identifier
        i_spectra <- lapply(seq_along(spectra), function(x) {
            if (inds_ev[[x]] != "-2") {
                sp <- spectra[[x]][spectra[[x]][, "PeakID"] == inds_ev[[x]],]
                ## if there is no matching MS/MS spectrum recorded set to NULL
                if (nrow(sp) == 0) sp <- NULL
            } else sp <- NULL
            return(sp)
        })
        
        ## if not all entries are NULL
        if (!is.null(unlist(i_spectra))) {
        
            ## retrieve the column "MSMS.spectrum" from each entry in i_spectra,
            ## this column contains information on the fragments and the 
            ## corresponding intensities, strsplit the entries and write them 
            ## to a matrix
            i_msms <- lapply(i_spectra, function(x) {
                if (!is.null(x)) {
                    msms <- strsplit(x[, "MSMS.spectrum"], split=" ")[[1]]  
                    msms <- lapply(msms, function(y) strsplit(y, split=":"))
                    msms <- matrix(unlist(msms), ncol=2, byrow=TRUE)
                    mode(msms) <- "numeric"
                    return(msms)
                }
            })
        
            ## bin within each spectra of i_spectra_msms
            ## sum up corresponding intensities from fragments that are  binned
            ## together
            i_msms_bin <- lapply(i_msms, function(x) {
                if (!is.null(x)) {
                    binSpectra(x, tol=tol, fun="sum")   
                }
            })
        
            ## deconvolute final spectra, bin across spectra
            deconv <- deconvolute(i_msms_bin, tol=tol)
        
        } else deconv <- NULL
        
        return(deconv)
    })
    
    ## construct names for res
    names(res) <- paste(sample_name, aln[, "Alignment.ID"], 
        aln[, "Average.Rt.min."], aln[, "Average.Mz"], sep="_")
    
    ## truncate res that it only contains the entries that are not NULL
    inds_keep <- unlist(lapply(res, function(x) !is.null(x)))
    res <- res[inds_keep]
    
    return(res)
}


#' @name binSpectra
#' @title Bin a spectra
#' @description \code{binSpectra} combines fragments that are close to each 
#' other
#' @usage binSpectra(spectra, tol)
#' @param spectra matrix, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#' @param tol numeric, tolerance parameter
#' @param fun character, either "max" or "sum"
#' @details Bin fragments and sum the intensities of corresponding fragments or 
#' take the maximum (depending on \code{fun}). If \code{fun} is equal to 
#' "max" the m/z value of the fragment with the highest intensity is 
#' retained togheter with the highest intensity, if \code{fun} is equal to 
#' "sum" the median of m/z values in the bin are retained and intensities 
#' are added. 
#' @return matrix containing in the first column m/z information and in the 
#' second column intensity information
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples
#' binSpectra(spectra, tol=0.01, fun="sum")
binSpectra <- function(spectra, tol=0.01, fun=c("max", "sum")) {
    
    ## match arguments for fun
    fun <- match.arg(fun)
    
    frag_s <- spectra[,1]
    steps <- (max(frag_s) - min(frag_s)) / tol
    bins <- tapply(frag_s, cut(frag_s, steps), mean)
    bins <- bins[!is.na(bins)]
    ## find nearest ones and sum all intensities up 
    inds <- lapply(frag_s, FUN = function(x) which.min(abs(x - bins)))
    inds <- unlist(inds)
    
    ## iterate through duplicated fragments and combine them
    for (j in names(which(table(inds) != 1))) {
        inds_dup <- which(inds == j)
        spectra_dup <- spectra[inds_dup,]
        
        ## either use max or sum the intensities
        if (fun == "max") {
            spectra[inds_dup, 1] <- spectra_dup[which.max(spectra_dup[, 2]), 1]
            spectra[inds_dup, 2] <- max(spectra_dup[, 2])    
            ## set all except the ones with the highest intensity to NA
            spectra[inds_dup[-which.max(spectra_dup[, 2])], ] <- NA 
        } 
        if (fun == "sum") {
            spectra[inds_dup, 1] <- median(spectra_dup[, 1])
            spectra[inds_dup, 2] <- sum(spectra_dup[, 2])   
            ## set all except the first one to NA
            spectra[inds_dup[-1],] <- NA 
        }
        
    }
    ## remove NA values
    spectra <- spectra[!is.na(spectra[,1]), ]
    return(spectra)
}

#' @name deconvolute
#' @title Deconvolute a final spectrum
#' @description \code{deconvolute} takes as input a list of MS/MS spectra. It will 
#' bin the m/z values and take the maximum of the corresponding intensity value.
#' @usage deconvolute(spectra, tol=0.01)
#' @param spectra list of matrices, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#' @param tol numeric, tolerance parameter for binning
#' @details Takes a list of matrices that contain spectral information and 
#' bins the m/z values. If several fragments are located in the same bin, 
#' the fragment with the highest intensity is retained and others are removed
#' from the spectrum. 
#' @return matrix containing in the first column m/z information and in the 
#' second column intensity information
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples
#' sp <- list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50)
#' deconvolute(spectra=sp, tol=0.01)
deconvolute <- function(spectra=list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50), tol=0.01) {
    ## rbind spectra
    spectra <- do.call("rbind", spectra)
    ## order spectra
    spectra <- spectra[order(spectra[,1]),]
    
    ## bin: take the max from the same corresonding fragment
    spectra <- binSpectra(spectra, tol=tol, fun="max")
    
    ## return
    return(spectra)
}

## use other function from MetCirc (>v1.15.0) and MSnbase (>v2.11.4) and 
## igraph (>v1.2.4.1)for visualization 
library(MSnbase)
library(igraph)

################################## workflow ####################################
## Maize MS2: create a similarity network
setwd("/home/thomas/Projects/Shijuan_MS2_Maize/")

## 1) 
## load file that contains alignment information
aln <- read.table("PeakID_0_2019820915.txt", sep="\t", fill=TRUE, skip=4, header=TRUE, stringsAsFactors=FALSE, quote="\"")
## truncate file
cols_keep <- c("Alignment.ID", "Average.Rt.min.", "Average.Mz", "Metabolite.name", "Adduct.type",
               "mixMK.ddms2.CID30", "mixMK.ddms2.CID40", "mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", 
               "mixMK.ddms2.HCD50", "mixML.ddms2.CID30", "mixML.ddms2.CID40", "mixML.ddms2.HCD30",             
               "mixML.ddms2.HCD40", "mixML.ddms2.HCD50")
aln <- aln[, cols_keep]
## keep only these rows that have alignment information
inds_remove <- apply(aln[,cols_keep[-c(1:5)]], 1, function(x) all(x == "-2"))
aln <- aln[!inds_remove, ]
## remove lines that are contain any NA
aln <- aln[!apply(aln, 1, function(x) any(is.na(x))), ]

## 2) 
## load files that contain spectra, do this for each sample individually
## for CID only 30 and 40 eV available, for HCD 30, 40 and 50 eV available
## for kernel
## CID
MK_cid30 <- read.table("mixMK-ddms2-CID30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_cid40 <- read.table("mixMK-ddms2-CID40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
## HCD
MK_hcd30 <- read.table("mixMK-ddms2-HCD30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd40 <- read.table("mixMK-ddms2-HCD40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd50 <- read.table("mixMK-ddms2-HCD50.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for leaf
## CID
ML_cid30 <- read.table("mixML-ddms2-CID30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_cid40 <- read.table("mixML-ddms2-CID40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
## HCD
ML_hcd30 <- read.table("mixML-ddms2-HCD30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd40 <- read.table("mixML-ddms2-HCD40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd50 <- read.table("mixML-ddms2-HCD50.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for sweet maize
## CID
swM_cid30 <- read.table("QC-ddms-CID30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_cid40 <- read.table("QC-ddms-CID40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
## HCD
swM_hcd30 <- read.table("QC-ddms-HCD30.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd40 <- read.table("QC-ddms-HCD40.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd50 <- read.table("QC-ddms-HCD50.txt", sep="\t", fill=TRUE, header=TRUE, stringsAsFactors=FALSE, quote="\"")

## 3) 
## truncate files and keep only columns that are of interest
cols_keep <- c("PeakID", "Title", "RT.min.", "Precursor.m.z", "Height", "Area", "MetaboliteName", "AdductIon", "Isotope", "S.N", "MSMS.spectrum")

MK_cid30 <- MK_cid30[MK_cid30[, "MSMS.spectrum"] != "", cols_keep]
MK_cid40 <- MK_cid40[MK_cid40[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd30 <- MK_hcd30[MK_hcd30[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd40 <- MK_hcd40[MK_hcd40[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd50 <- MK_hcd50[MK_hcd50[, "MSMS.spectrum"] != "", cols_keep]
ML_cid30 <- ML_cid30[ML_cid30[, "MSMS.spectrum"] != "", cols_keep]
ML_cid40 <- ML_cid40[ML_cid40[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd30 <- ML_hcd30[ML_hcd30[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd40 <- ML_hcd40[ML_hcd40[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd50 <- ML_hcd50[ML_hcd50[, "MSMS.spectrum"] != "", cols_keep]

## 4) 
## truncate files that they only contain the spectra that are aligned
## (spectra that are present in aln)
MK_cid30 <- MK_cid30[MK_cid30[, "PeakID"] %in% aln[, "mixMK.ddms2.CID30"],]
MK_cid40 <- MK_cid40[MK_cid40[, "PeakID"] %in% aln[, "mixMK.ddms2.CID40"],]

MK_hcd30 <- MK_hcd30[MK_hcd30[, "PeakID"] %in% aln[, "mixMK.ddm2.HCD30"],]
MK_hcd40 <- MK_hcd40[MK_hcd40[, "PeakID"] %in% aln[, "mixMK.ddms2.HCD40"],]
MK_hcd50 <- MK_hcd50[MK_hcd50[, "PeakID"] %in% aln[, "mixMK.ddms2.HCD50"],]

ML_cid30 <- ML_cid30[ML_cid30[, "PeakID"] %in% aln[, "mixML.ddms2.CID30"],]
ML_cid40 <- ML_cid40[ML_cid40[, "PeakID"] %in% aln[, "mixML.ddms2.CID40"],]

ML_hcd30 <- ML_hcd30[ML_hcd30[, "PeakID"] %in% aln[, "mixML.ddm2.HCD30"],]
ML_hcd40 <- ML_hcd40[ML_hcd40[, "PeakID"] %in% aln[, "mixML.ddms2.HCD40"],]
ML_hcd50 <- ML_hcd50[ML_hcd50[, "PeakID"] %in% aln[, "mixML.ddms2.HCD50"],]


## 5) 
## assemble spectra
assembly_MK_hcd <- assemblySpectra(spectra=list(MK_hcd30, MK_hcd40, MK_hcd50), 
    aln=aln, cols_aln=c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD")

## 6) 
## convert final deconvoluted spectra to Spectrum2 objects from MSnbase

## load the MSnbase package
library(MSnbase)

## retrieve retention time, precursorMZ, m/z of fragments and corresponding 
## intensities from the assembled spectrum
rt_MK_hcd <- unlist(lapply(strsplit(names(assembly_MK_hcd), split="_"), "[", 4))
rt_MK_hcd <- as.numeric(rt_MK_hcd)
prec_mz_MK_hcd <- unlist(lapply(strsplit(names(assembly_MK_hcd), split="_"), "[", 5))
prec_mz_MK_hcd <- as.numeric(prec_mz_MK_hcd)
mz_MK_hcd <- lapply(assembly_MK_hcd, function(x) x[,1])
int_MK_hcd <- lapply(assembly_MK_hcd, function(x) x[,2])

## create Spectrum2 objects with the corresponding information
spN_l <- lapply(1:length(assembly_MK_hcd), function(x) new("Spectrum2", 
            rt=rt_MK_hcd[x], precursorMz=prec_mz_MK_hcd[x], 
            mz=mz_MK_hcd[[x]], intensity=int_MK_hcd[[x]]))

## 7)
## create Spectra object from the list of Spectrum2 objects
spectra_li <- MSnbase::Spectra(spN_l, elementMetadata=S4Vectors::DataFrame(show=rep(TRUE, length(spN_l)))) 

## 8) 
## calculate pair-wise similarities using the normalizeddotproduct
similarityMat <- compare_Spectra(spectra_li, fun=normalizeddotproduct, binSize=0.01)

## 9) 
## visualize network 
library(igraph)
colnames(similarityMat) <- rownames(similarityMat) <- names(assembly_MK_hcd)
graph_from_adjacency_matrix(similarityMat, weighted = TRUE, mode = "undirected")
