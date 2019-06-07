## helper functions for normalization.R
## taken from https://bitbucket.org/veitveit/vsclust/src/master/
##library(e1071FuzzVec)

statWrapper <- function(dat, NumReps, NumCond, isPaired=F, isStat) {
    
    # extend to 702 cases:
    LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
    
    qvals <- statFileOut <- Sds <- NULL
    fdat <- dat
    tdat <- NULL
    print(isStat)
    if (isStat) {
        if(ncol(dat)!=NumReps*NumCond)
            print("Number of data columns must correspond to product of conditions and replicates!")
        if (isPaired) {
            print("Paired")
            ttt <- SignAnalPaired(dat, NumCond, NumReps)
        } else {
            ttt <- SignAnal(dat, NumCond, NumReps)
        }
        
        tdat <- rowMeans(dat[, seq(1, NumReps*NumCond, NumCond)], na.rm=T) ## was 1, reps*cond, cond
        
        print(Sds)
        
        Sds <- ttt$Sds
        qvals <- ttt$plvalues
        colnames(qvals) <- paste("qvalue ", LETTERS702[2:(NumCond)], "vsA", sep="")
        for (i in 2:NumCond) {
            tdat <- cbind(tdat, rowMeans(dat[, seq(i, NumReps*NumCond, NumCond)], na.rm=T))
        }
        colnames(tdat)<-paste("Mean of log ", LETTERS702[1:(NumCond)], sep="")
        dat <- cbind(tdat, Sds=Sds)
        #      dat <- dat[rowSums(is.na(ttt$plvalues))==0,]
        tdat <- dat[, 1:(ncol(dat)-1)]
        Sds <- dat[, ncol(dat)]
    } else {
        Sds <- dat[, ncol(dat)]
        tdat <- dat[, 1:(ncol(dat)-1)]
        NumReps <- 1
        NumCond <- ncol(dat) - 1
    }
    dat[, ncol(dat)] <- dat[, ncol(dat)] / rowSds(as.matrix(tdat), na.rm=T)
    
    if (isStat) {
        statFileOut <- cbind(dat[, 1:(ncol(dat)-1)], Sds, qvals)
    } else {
        statFileOut <- cbind(dat[, 1:(ncol(dat)-1)], Sds)
    }
    
    # Remove columns with only 5% of the data
    plab <- rep(1:NumCond,NumReps)
    plab <- plab[colSums(!is.na(fdat)) > 0.05*nrow(fdat)]
    pcaDat <- fdat[, colSums(!is.na(fdat)) > 0.05*nrow(fdat)]
    
    # pcaDat <- t(pcaDat[complete.cases(pcaDat),])
    pcaDat <- (pcaDat[complete.cases(pcaDat), -c(1:3)])
    print(which(apply(pcaDat, 2, sd) == 0))
    ##validate(need(length(pcaDat)> 0, "Principal component analysis not shown as too many missing values"))      
    ##validate(need(nrow(pcaDat)> 10, "Principal component analysis not shown as too many missing values"))      
    pca <- prcomp(pcaDat, scale=T, retx=T)
    # scores <- pca$x
    # loadings <- pca$rotation
    scores <- pca$rotation
    loadings <- pca$x
    par(mfrow=c(1, 3))
    tSds <- ifelse(is.na(dat[, ncol(dat)]), max(dat[, ncol(dat)], na.rm=T), sqrt(dat[, ncol(dat)]))
    print(max(as.integer(255/max(1/tSds)/tSds)))
    plot(loadings, cex=tSds, pch=16, 
         col=paste("#000000",sprintf("%02X", as.integer(255/max(1/tSds)/tSds)), sep=""))
    title(main="Principal component analysis of data set (loadings)", sub="The point size corresponds to the estimated standard deviation")
    plot(scores, pch=19,col=rainbow(NumCond)[plab])
    title(main="Principal component analysis of data set (scores)",sub="Colors denote different conditions")
    plot.new()
    legend("topright", paste("Condition", 1:NumCond), pch=rep(19, NumCond), col=rainbow(NumCond)[1:NumCond])
    
    ## Preparing output
    Out <- list(dat=dat, qvals=qvals, statFileOut=statFileOut)
    Out
    
}

SignAnal <- function(Data,NumCond,NumReps) {
    ##########################################################
    # significance analysis
    Reps <- rep(1:NumCond,NumReps)
    design <- model.matrix(~0+factor(Reps-1))
    colnames(design) <- paste("i", c(1:NumCond),sep="")
    contrasts<-NULL
    First <- 1
    for (i in (1:NumCond)[-First]) contrasts <- append(contrasts, paste(colnames(design)[i], "-", colnames(design)[First], sep=""))
    contrast.matrix <- makeContrasts(contrasts=contrasts, levels=design)
    print(dim(Data))
    lm.fitted <- lmFit(Data, design)
    
    lm.contr <- contrasts.fit(lm.fitted, contrast.matrix)
    lm.bayes<-eBayes(lm.contr)
    topTable(lm.bayes)
    plvalues <- lm.bayes$p.value
    qvalues <- matrix(NA, nrow=nrow(plvalues), ncol=ncol(plvalues), dimnames=dimnames(plvalues))
    # qvalue correction
    for (i in 1:ncol(plvalues)) {
        tqs <- qvalue(na.omit(plvalues[, i]))$qvalues
        qvalues[names(tqs), i] <- tqs
    }
    return(list(plvalues=qvalues, Sds=sqrt(lm.bayes$s2.post)))
}

SignAnalPaired <- function(Data, NumCond, NumReps) {
    ##########################################################
    # significance analysis
    MAData<-Data[, 2:(NumCond)]-Data[, 1]
    for (i in 1:(NumReps-1))
        MAData <- cbind(MAData, Data[, (i*NumCond+1)+1:(NumCond-1)]-Data[, (i*NumCond+1)])
    rownames(MAData) <- rownames(Data)
    MAReps <- rep(1:(NumCond-1), NumReps)
    ##limma with ratios
    design <- plvalues <- NULL
    for (c in (1:(NumCond-1))) {
        design <- cbind(design, as.numeric(MAReps==c))
    }
    lm.fittedMA <- lmFit(MAData, design)
    lm.bayesMA <- eBayes(lm.fittedMA)
    topTable(lm.bayesMA)
    plvalues <- lm.bayesMA$p.value
    qvalues <- matrix(NA, nrow=nrow(plvalues), ncol=ncol(plvalues), dimnames=dimnames(plvalues))
    # qvalue correction
    for (i in 1:ncol(plvalues)) {
        tqs <- qvalue(na.omit(plvalues[, i]))$qvalues
        qvalues[names(tqs), i] <- tqs
    }
    
    return(list(plvalues=qvalues, Sds=sqrt(lm.bayesMA$s2.post)))
}

ClustComp <- function(tData, NSs=100, NClust=NClust, Sds=Sds, cores=1) {
    D <- ncol(tData)
    d <- sqrt(D/2)
    dims <- dim(tData)                                                                                                           
    mm <- 1 + (1418/dims[1] + 22.05)* dims[2]^(-2) + (12.33/dims[1] + 0.243)*dims[2]^(-0.0406*log(dims[1])-0.1134)
    
    ### d_i and d_t
    difunc <- function(c,D) { x <- 0:c; sum(choose(c, x) / (x*D+1)*(-1)^x) }
    
    di <- difunc(NClust, D)  / sqrt(pi) * gamma(D/2+1)^(1/D)
    dt <- NClust^(-1/D)
    
    p <- dnorm(di, 0, Sds) * (1-dnorm(dt, 0, Sds))^(NClust - 1)
    
    m <- mm + p*mm*(D/3-1)
    m[m==Inf] <- 0
    m[m==0] <- NA
    m[is.na(m)] <- mm*10
    #m<-rowMaxs(cbind(m,mm,na.rm=T))
    
    ## If m for highest Sd is mm then all = mm
    if (m[which.max(Sds)]== mm) 
        m[1:length(m)] <- mm
    
    #   plot(Sds,m,cex=0.5,pch=15,col=rainbow(maxClust)[NClust])
    #   abline(h=mm)
    
    colnames(tData)<-NULL
    PExpr <- new("ExpressionSet", expr=as.matrix(tData)) 
    PExpr.r <- filter.NA(PExpr, thres=0.25)
    PExpr <- fill.NA(PExpr.r, mode = "mean")
    tmp <- filter.std(PExpr, min.std=0, visu=F)
    PExpr2 <- standardise(PExpr)
    
    ## VSclust
    cls <- mclapply(1:NSs, function(x) cmeans(exprs(PExpr2), NClust, m=m, verbose=F, iter.max=1000), mc.cores=cores)
    print(cls[[1]])
    Bestcl <- cls[[which.min(lapply(cls, function(x) x$withinerror))]]
    ## Standard clust
    cls <- mclapply(1:NSs, function(x) cmeans(exprs(PExpr2), NClust, m=mm, verbose=F, iter.max=1000), mc.cores=cores)
    Bestcl2 <- cls[[which.min(lapply(cls, function(x) x$withinerror))]]
    
    # return validation indices
    list(indices=c(min(dist(Bestcl$centers)),  cvalidate.xiebeni(Bestcl, mm), ## VSclust
                   min(dist(Bestcl2$centers)), cvalidate.xiebeni(Bestcl2, mm)), ## standard Clust
         Bestcl=Bestcl, Bestcl2=Bestcl2, m=m, withinerror=Bestcl$withinerror, withinerror2=Bestcl2$withinerror) 
}

## Wrapper for estimation of cluster number
estimClustNum <- function(dat, maxClust=25, cores=1) {
    # print(head(rowSds(as.matrix(dat[,1:(ncol(dat)-1)]),na.rm=T)))
    ##Sds <- dat[,ncol(dat)] / rowSds(as.matrix(dat[,1:(ncol(dat)-1)]),na.rm=T)
    ClustInd <- matrix(NA, nrow=maxClust, ncol=6)
    multiOut <- mclapply(3:maxClust, function(x) {    
        if (!is.null(getDefaultReactiveDomain())) {
            incProgress(1, detail = paste("Running cluster number",x))
        } else {
            print(paste("Running cluster number",x))
        }
        clustout <- ClustComp(dat[,1:(ncol(dat)-1)], NClust=x, Sds=dat[,ncol(dat)], NSs=100, cores)
        c(clustout$indices, sum(rowMaxs(clustout$Bestcl$membership)>0.5), 
          sum(rowMaxs(clustout$Bestcl2$membership)>0.5))
    }, mc.cores=cores)
    for (NClust in 3:maxClust) 
        ClustInd[NClust, ] <- multiOut[[NClust-2]]

    dmindist <- c(which.max(ClustInd[3:(maxClust-2), 1]-ClustInd[4:(maxClust-1), 1])+2,
                  which.max(ClustInd[3:(maxClust-2), 3]-ClustInd[4:(maxClust-1), 3])+2)
    dxiebeni <- c(which.min(ClustInd[3:(maxClust-1), 2])+2,
                  which.min(ClustInd[3:(maxClust-1), 4])+2)  
    print(dmindist)
    plot.new() ## clean up device
    par(mfrow=c(1,3))
    plot(3:maxClust, ClustInd[3:(maxClust), 1],col=2 , type="b", 
         main="Min. centroid distance\n(Highest jump is best)", xlab="Number of clusters", 
         ylab="Index",ylim=c(min(ClustInd[, c(1, 3)], na.rm=T), max(ClustInd[, c(1, 3)],na.rm=T)))
    lines(3:maxClust, ClustInd[3:(maxClust), 3],col=3, type="b")
    points(dmindist[1], ClustInd[dmindist[1], 1], pch=15, col=1, cex=2)
    legend("topright", legend=c("VSClust","Standard"), lty=c(1, 1), col=2:3)
    grid(NULL, NA, lwd=1, col=1)
    plot(3:maxClust, ClustInd[3:(maxClust), 2], col=2, type="b", 
         main="Xie-Beni index\n(Lowest is best)", xlab="Number of clusters", 
         ylab="Index",ylim=c(min(ClustInd[, c(2, 4)], na.rm=T), max(ClustInd[, c(2,4)], na.rm=T)))
    lines(3:maxClust, ClustInd[3:(maxClust), 4],type="b", col=3)
    points(dxiebeni[1], ClustInd[dxiebeni[1], 2], pch=15, col=1, cex=2)
    legend("topright",legend=c("VSClust","Standard"), lty=c(1, 1), col=2:3)
    grid(NULL, NA, lwd=1, col=1)
    plot(3:maxClust, ClustInd[3:(maxClust), 5], col=2, type="b", 
         main="Total number of assigned features", xlab="Number of clusters", 
         ylab="Assigned features", ylim=c(min(ClustInd[, 5:6], na.rm=T), max(ClustInd[, 5:6], na.rm=T)))
    lines(3:maxClust, ClustInd[3:(maxClust), 6], type="b", col=3)
    legend("topright",legend = c("VSClust","Standard"), lty=c(1,1),col=2:3)
    # finally plot
    p <- recordPlot()
    
    # Output
    Out <- list(ClustdInd=ClustInd, p=p, 
                numclust_vs_dmindist=dmindist[1], numclust_st_dmindist=dmindist[2], 
                numclust_vs_dxiebeni=dxiebeni[1], numclust_st_dxiebeni=dxiebeni[2])
    Out
}

SwitchOrder <- function(Bestcl, NClust) {
    switching <- table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5])
    switching <- as.numeric(names(sort(switching, decreasing=T)))
    if (length(switching) < NClust) 
        switching <- c(switching, which(!((1:NClust) %in% switching)))
    switching2 <- 1:NClust
    names(switching2) <- switching
    tBest <- Bestcl
    tBest$centers <- Bestcl$centers[switching, ]
    rownames(tBest$centers) <- 1:NClust
    tBest$size <- Bestcl$size[switching]
    tBest$cluster <- switching2[as.character(Bestcl$cluster)]
    names(tBest$cluster) <- names(Bestcl$cluster)
    tBest$membership <- Bestcl$membership[, switching]
    colnames(tBest$membership) <- 1:NClust
    tBest
}

runClustWrapper <- function(dat, NClust, proteins=NULL, VSClust=T, cores) {
    # dat <- dat[rowSums(is.na(dat))==0,]
    PExpr <- new("ExpressionSet", expr=as.matrix(dat[, 1:(ncol(dat)-1)])) ## was ncol(dat)-1
    PExpr.r <- filter.NA(PExpr, thres=0.25)
    PExpr <- fill.NA(PExpr.r, mode = "mean")
    tmp <- filter.std(PExpr, min.std=0, visu=F)
    PExpr <- standardise(PExpr)
    
    
    clustout <- ClustComp(exprs(PExpr), NClust=NClust, Sds=dat[, ncol(dat)], NSs=100, cores)
    if (VSClust) {
        Bestcl <- clustout$Bestcl
    } else {
        Bestcl <- clustout$Bestcl2
    }
    Bestcl <- SwitchOrder(Bestcl, NClust)
    ##########

    # sorting for membership values (globally)
    Bestcl$cluster <- Bestcl$cluster[order(rowMaxs(Bestcl$membership, na.rm=T))]
    Bestcl$membership <- Bestcl$membership[order(rowMaxs(Bestcl$membership, na.rm=T)), ]
    PExpr <- PExpr[names(Bestcl$cluster), ]
    
    if (!is.null(getDefaultReactiveDomain()))
        incProgress(0.7, detail=paste("Plotting", NClust))
    
    plot.new() ## clean up device
    par(lwd=0.25)
    oldmar <- par("mar")
    par(mar=c(2,2,3,3),mgp=c(2, 1, 0))
    par(mar=par("mar") / max(1, NClust/20))
    mfuzz.plot2(PExpr, cl=Bestcl, mfrow=c(round(sqrt(NClust)), ceiling(sqrt(NClust))), min.mem=0.5, x11=F, colo="fancy")
    # Mfuzz::mfuzzColorBar(col="fancy")
    p <- recordPlot()
    par(lwd=1,mar=oldmar)
    
    colnames(Bestcl$membership) <- paste("membership of cluster", colnames(Bestcl$membership))
    outFileClust <- exprs(PExpr)
    if (!is.null(proteins)) {
        outFileClust <- cbind(outFileClust, names=as.character(proteins[rownames(outFileClust)]))
    }
    
    rownames(Bestcl$centers) <- paste("Cluster",rownames(Bestcl$centers))
    ClustInd <- as.data.frame(table(Bestcl$cluster[rowMaxs(Bestcl$membership)>0.5]))
    if (ncol(ClustInd) == 2)
        colnames(ClustInd) <- c("Cluster","Members")
    else 
        ClustInd <- cbind(1:max(Bestcl$cluster), rep(0,max(Bestcl$cluster)))
    
    ## Output
    Out <- list(dat=PExpr, Bestcl=Bestcl, p=p, outFileClust=outFileClust, ClustInd=ClustInd)
    return(Out)
}

cvalidate.xiebeni <- function(clres,m) {                                                
    xrows <- dim(clres$me)[1]                                                
    minimum <- -1                                                            
    error <- clres$within                                                    
    ncenters <- dim(clres$centers)[1]                                        
    for (i in 1:(ncenters - 1)) {                                            
        for (j in (i + 1):ncenters) {                                        
            diff <- clres$ce[i, ] - clres$ce[j, ]                            
            diffdist <- t(diff) %*% t(t(diff))                               
            if (minimum == -1)                                               
                minimum <- diffdist                                            
            if (diffdist < minimum)                                          
                minimum <- diffdist                                            
        }                                                                    
    }                                                                        
    xiebeni <- error/(xrows * minimum)                                       
    return(xiebeni)                                                          
}

