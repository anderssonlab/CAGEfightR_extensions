require("CAGEfightR")

## Object: GRanges
## ctss: RangedSummarizedExperiment
CTSSvariation <- function(object, ctss, inputAssay="counts", outputColumn="MSE", fn=mean_squared_error) {

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate variation
    message("Calculating variation...")
    
    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids, function(x) {
        fn(data[x,])
    }))

    mcols(object)[, outputColumn] <- value

    object
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

random_mean_squared_error <- function(d) {
    d <- as.matrix(d)
    t <- rowSums(d)
    n <- nrow(d)
    p <- colSums(d) / sum(t)
    
    mean(sapply(1:10, function(i) {
        ## Subsample the pooled counts to each sample's proportion of the total
        d <- scale(sapply(p, function(x) rbinom(n, t, x)), center=TRUE, scale=TRUE)
        cent <- rowMeans(d,na.rm=TRUE)
        mean(apply(d, 2, function(x) sum((x - cent)^2, na.rm=TRUE)) / nrow(d), na.rm=TRUE)
    }), na.rm=TRUE)
}

random_weighted_mean_squared_error <- function(d) {
    d <- as.matrix(d)
    t <- rowSums(d)
    n <- nrow(d)
    p <- colSums(d) / sum(t)
    
    mean(sapply(1:10, function(i) {
        ## Subsample the pooled counts to each sample's proportion of the total
        d <- scale(sapply(p, function(x) rbinom(n, t, x)), center=TRUE, scale=TRUE)
        cent <- rowMeans(d,na.rm=TRUE)
        sum(p*apply(d, 2, function(x) sum((x - cent)^2, na.rm=TRUE)) / nrow(d), na.rm=TRUE)
    }), na.rm=TRUE)
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

random_divergence <- function(d) {
    d <- as.matrix(d)
    t <- rowSums(d)
    n <- nrow(d)
    p <- colSums(d) / sum(t)

    mean(sapply(1:10, function(i) {
        ## Subsample the pooled counts to each sample's proportion of the total
        d <- sapply(p, function(x) rbinom(n, t, x))
        
        t <- rowMeans(d)
        cent <- t / sum(t)
        
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
    }), na.rm=TRUE)
}
