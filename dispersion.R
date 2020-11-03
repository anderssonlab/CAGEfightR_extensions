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

madm <- function(x) {
    m <- median(x,na.rm=TRUE)
    if (m==0)
        return(NA)
    median(abs(x-m),na.rm=TRUE)/m
}

tss_madm <- function(d, pseudo=NULL, max_zero_prop=0.5) {
    d <- as.matrix(d)
    if (max_zero_prop < 1)
    {
        propz <- apply(d,1,function(x) sum(x==0))/ncol(m)
        idx <- which(propz <= max_zero_prop)
    }
    else
        idx <- 1:nrow(d)
    if (!is.null(pseudo))
        d <- t(t(d)+pseudo)
    meds <- apply(d,1,median)
    idx <- idx[which(meds[idx]>0)]
    dnorm <- apply(d,2,function(x) x/sum(x))
    median(apply(dnorm[idx,,drop=FALSE], 1, madm),na.rm=TRUE)
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
