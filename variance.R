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
    value <- unlist(bplapply(ids, function(x) {
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

