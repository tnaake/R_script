#' @name cellColor
#' 
#' @title cellColor - helper function for metIDQ_get_high_quality_features
cellColor <- function(style) {
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) "")
    rgb <- paste(rgb, collapse = "")
    return(rgb)
}


#' @name metIDQ_getHighQualityFeatures
#'
#' @title Get a vector with the high-quality features to retain
metIDQ_getHighQualityFeatures <- function(file, threshold = 0.66) {
    ## get the background color
    wb <- xlsx::loadWorkbook(file = file)
    sheet1 <- xlsx::getSheets(wb)[[1]]
    rows <- xlsx::getRows(sheet1)
    cells <- xlsx::getCells(rows)
    styles <- sapply(cells, xlsx::getCellStyle)
    bg <- sapply(styles, cellColor) 
    
    ## use values to get the indices where the actual data is stored
    values <- sapply(cells, xlsx::getCellValue)
    
    ## convert bg and values to matrix that it has the same dimension as before
    if (all(names(bg) != names(values))) stop("dimensions do not match")
    row_index <- as.numeric(unlist(lapply(strsplit(names(bg), split = "[.]"), "[", 1)))
    col_index <- as.numeric(unlist(lapply(strsplit(names(bg), split = "[.]"), "[", 2)))
    mat_values <- mat_bg <- matrix("", ncol = max(col_index), nrow = max(row_index))
    for (i in 1:length(bg)) {
        mat_bg[row_index[i], col_index[i]] <- bg[i]
        mat_values[row_index[i], col_index[i]] <- values[i]
    }
    
    ## obtain the indices of the cells where the values/background colors are 
    ## stored in
    row_inds <- which(!is.na(mat_values[, 1]))
    row_inds <- row_inds[!row_inds %in% 1:2]
    col_inds <- rbind(which(mat_values == "C0", arr.ind = TRUE),
                      which(mat_values == "Choline", arr.ind = TRUE))
    colnames(mat_values) <- colnames(mat_bg) <- mat_values[col_inds[1, "row"], ]
    col_inds <- seq(col_inds[1, "col"], col_inds[2, "col"])
    
    ## iterate through the columns of mat_bg and check the color values 
    ## against the QC of MetIDS## "00cd66" == green, "87ceeb" == lightblue
    valid <- apply(mat_bg[row_inds, col_inds], 2, 
                   function(x) sum(x %in% c("00cd66", "87ceeb"))) 
    
    ## require that at least threshold*100% values per metabolite are 
    ## "green"/"lightblue", 
    valid / length(row_inds)  > threshold
}