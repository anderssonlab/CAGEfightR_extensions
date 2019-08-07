require("CAGEfightR")

## Object: GRanges
## ctss: RangedSummarizedExperiment
CTSSvariation <- function(object, ctss, inputAssay="counts", outputColumn="MSE", standardize=TRUE, fn=mean_squared_error) {

    ## Find overlaps
    message("Finding overlaps...")

    hits <- findOverlaps(query = ctss,
                         subject = object,
                         select = "arbitrary")
    hits <- factor(hits, levels=seq_along(object))
    ids <- as.numeric(unlist(split(1:length(hits), hits)))
    ctss <- ctss[ids,]
    ids <- split(1:nrow(ctss), hits[-which(is.na(hits))])

    ## Calculate variances
    message("Calculating variation...")
    
    data <- assay(ctss, inputAssay)
    value <- unlist(bplapply(ids[1:11], function(x) {
        d <- scale(data[x,], center=standardize, scale=standardize)
        fn(d)
    }))

    mcols(object)[, outputColumn] <- value

    object
}

mean_squared_error <- function(d) {
    val <- which(apply(d, 2, function(x) !(any(is.na(x)))))
    cent <- rowMeans(d[,val])
    sum(apply(d[,val], 2, function(x) sum((x - cent)^2))) / nrow(d)
}

average_divergence <- function(d) {
    val <- which(apply(d, 2, function(x) !(any(is.na(x)))))
    cent <- rowMeans(d[,val])
    cent <- cent / sum(cent)
    
    sum(apply(d[,val], 2, function(x) {
        x <- x / sum(x)
        m <- 0.5 * (cent + x)
        z <- 0.5 * (x * log2(x/m) + cent * log2(cent/m))
        
        z[z == -Inf] <- 0
        z[z == Inf] <- 0
        z[is.na(z)] <- 0
        sqrt(sum(z))
    })) / length(val)
}

