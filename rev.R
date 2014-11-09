

reverse <- function(x) {
    x <- strsplit(x, split = NULL)
    x <- x[[1]]
    y <- x
    y[which(x == "A")] <- "T"
    y[which(x == "T")] <- "A"
    y[which(x == "G")] <- "C"
    y[which(x == "C")] <- "G"
    y[length(y):1] <- y[1:length(y)]
    y <- paste0(y, collapse="")
    return(y)
}
