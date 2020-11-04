require("CAGEfightR")

## Object: GRanges
## ctss: RangedSummarizedExperiment
CTSSdispersion <- function(object, ctss, inputAssay="counts", outputColumn="TSS_MADM", fn=tss_madm, ...) {

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate dispersion
    message("Calculating TSS usage dispersion...")

    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids, function(x) {
        fn(data[x,],...)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
TCdispersion <- function(object, inputAssay="counts", outputColumn="TC_MADM", fn=tc_madm, ...) {

    message("Calculating TC dispersion...")

    data <- assay(object, inputAssay)
    value <- unlist(bplapply(1:nrow(data), function(i) {
        fn(data[i,],...)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
calcMedian <- function(object, inputAssay="counts", outputColumn="median") {

    data <- assay(object, inputAssay)
    value <- apply(data,1,median,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

## object: RangedSummarizedExperiment
calcMean <- function(object, inputAssay="counts", outputColumn="mean") {

    data <- assay(object, inputAssay)
    value <- rowMeans(data,na.rm=TRUE)

    mcols(object)[, outputColumn] <- value

    object
}

## Object: GRanges
## ctss: RangedSummarizedExperiment
numCTSSs <- function(object, ctss, inputAssay="counts", outputColumn="num_CTSSs", max_zero_prop=0.5) {

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate dispersion
    message("Calculating number of CTSSs with valid zero proportion...")

    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids, function(x) {
        x <- as.matrix(x)
        sum(apply(x,1,function(y) sum(y==0))/ncol(x) <= max_zero_prop)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

madm <- function(x) {
    m <- median(x,na.rm=TRUE)
    if (m==0)
        return(NA)
    median(abs(x-m),na.rm=TRUE)/m
}

tss_madm <- function(d, pseudo=rep(0,ncol(d)), max_zero_prop=0.5) {
    d <- as.matrix(d)
    if (nrow(d)==1)
        return(0)
    if (max_zero_prop < 1)
    {
        propz <- apply(d,1,function(x) sum(x==0))/ncol(d)
        idx <- which(propz <= max_zero_prop)
    }
    else
        idx <- 1:nrow(d)
    if (length(idx)==0)
        return(NA)
    if (length(idx)==1)
        return(0)
    d <- t(t(d)+pseudo)
    dnorm <- matrix(0,nrow=nrow(d),ncol=ncol(d))
    dnorm[idx,] <- apply(d[idx,,drop=FALSE],2,function(x) x/sum(x))
    median(apply(dnorm[idx,,drop=FALSE], 1, madm),na.rm=TRUE)
}

tc_madm <- function(d, pseudo=rep(0,length(d)), max_zero_prop=0.5) {
    d <- as.numeric(d)
    if (max_zero_prop < 1)
    {
        propz <- sum(d==0)/length(d)
        if (propz > max_zero_prop)
            return(NA)
    }
    d <- d+pseudo
    madm(d)
}

mean_squared_error <- function(d) {
    d <- as.matrix(d)
    d <- scale(d, center=TRUE, scale=TRUE)
    cent <- rowMeans(d,na.rm=TRUE)
    mean(apply(d, 2, function(x) sum((x - cent)^2)) / nrow(d), na.rm=TRUE)
}

weighted_mean_squared_error <- function(d) {
    d <- as.matrix(d)
    p <- colSums(d) / sum(d)
    d <- scale(d, center=TRUE, scale=TRUE)
    cent <- rowMeans(d,na.rm=TRUE)
    sum(p*apply(d, 2, function(x) sum((x - cent)^2)) / nrow(d), na.rm=TRUE)
}

average_divergence <- function(d) {
    d <- as.matrix(d)
    cent <- rowMeans(d)
    cent <- cent / sum(cent)

    ## Average Jensen-Shannon divergence versus the cluster centroid
    mean(apply(d, 2, function(x) {
        x <- x / sum(x)
        m <- 0.5 * (cent + x)
        z <- 0.5 * (x * log2(x/m) + cent * log2(cent/m))

        z[z == -Inf] <- 0
        z[z == Inf] <- 0
        z[is.na(z)] <- 0
        sqrt(sum(z))
    }))
}
