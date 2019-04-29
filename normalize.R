require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## object: SummarizedExperiment
## inputAssay: assay to normalize
## outputAssay: where to store normalized results
## conditionalColumn: numeric row-wise vector to normalize against
## offsetAssay: where to store offset

## GC normalization based on approach described in Pickrell et al 2010, Nature

conditionalNormalize <- function(object, inputAssay="counts", outputAssay="normalized", conditionalColumn="GC", offsetAssay=NULL, bins=200, sizeFactors=NULL, aggregate.fn=sum) {
    
    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                conditionalColumn %in% colnames(mcols(object)))

    message("binning data according to conditional column...")
    
    y <- assay(object,inputAssay)
    x <- as.numeric(mcols(object)[,conditionalColumn])
    b <- as.character(Hmisc::cut2(x, unique(quantile(x, seq(1/bins, 1, 1/bins)))))
    b.m <- aggregate(x, by=list(b), FUN=mean)
    n <- b.m[,1]
    b.m <- b.m[,2]
    names(b.m) <- n

    y.b <- aggregate(as.matrix(y),by=list(b),FUN=aggregate.fn)
    b.g <- y.b[,1]
    rownames(y.b) <- b.g
    y.b <- y.b[names(b.m),]
    y.b <- data.matrix(y.b[,-1])

    message("calculating relative enrichments...")

    if (is.null(sizeFactors))
        sizeFactors <- rep(1,ncol(y))

    y.b <- sapply(1:ncol(y.b), function(i) y.b[,i] / sizeFactors[i])
    s <- colSums(y.b) / sum(y.b)
    
    f <- log2( t( t(y.b / rowSums(y.b)) / s ) )
    fit <- lapply(1:ncol(y), function(i) {
        b.m.i <- b.m
        f.i <- f[,i]
        missing <- which(is.na(f[,i]) | is.infinite(f[,i]))
        if (length(missing)>0) {
            b.m.i <- b.m.i[-missing]
            f.i <- f.i[-missing]
        }
        smooth.spline(b.m.i, f.i, spar=1)
    })

    message("normalizing data according to enrichments...")

    offset <- Matrix::Matrix(sapply(1:ncol(y), function(i) predict(fit[[i]],b.m[b])$y))
    y_normed <- Matrix::Matrix(ceiling(y * 2^(-offset)))

    dimnames(y_normed) <- dimnames(y)
    dimnames(offset) <- dimnames(y)
    
    assay(object, outputAssay) <- y_normed
    if (!is.null(offsetAssay))
        assay(object, offsetAssay) <- offset

    ## return
    object
}

normalizeBySizeFactors <- function(object, sizeFactors, inputAssay="counts", outputAssay="normalized") {

    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                length(sizeFactors) == ncol(object))

    assay(object, outputAssay) <- Matrix::t(Matrix::t(assay(object, inputAssay)) / sizeFactors)

    object
}
