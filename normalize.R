require("CAGEfightR")
require("assertthat")

source("CAGEfightR_extensions/utils.R")

## object: SummarizedExperiment
## inputAssay: assay to normalize
## outputAssay: where to store normalized results
## conditionalColumn: numeric row-wise vector to normalize against

## GC normalization approach described in Pickrell et al 2010, Nature

conditionalNormalize <- function(object, inputAssay="counts", outputAssay="normalized", conditionalColumn="GC", bins=200) {
    
    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                conditionalColumn %in% colnames(mcols(object)))

    message("binning data according to conditional column...")
    
    y <- assay(object,inputAssay)
    x <- as.numeric(mcols(object)[,conditionalColumn])
    b <- as.character(Hmisc::cut2(x,g=bins,labels=FALSE))
    b.m <- aggregate(x, by=list(b), FUN=mean)
    n <- b.m[,1]
    b.m <- b.m[,2]
    names(b.m) <- n

    y.b <- aggregate(as.matrix(y),by=list(b),FUN=sum)
    b.g <- y.b[,1]
    rownames(y.b) <- b.g
    y.b <- y.b[names(b.m),]
    y.b <- data.matrix(y.b[,-1])

    message("calculating relative enrichments...")
    
    f <- log2( t( t(y.b / rowSums(y.b)) / (colSums(y.b) / sum(y.b)) ) )
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
    
    y_normed <- Matrix::Matrix(sapply(1:ncol(y), function(i) ceiling(y[,i] * 2^(-1*predict(fit[[i]],b.m[b])$y))))

    dimnames(y_normed) <- dimnames(y)
    
    assay(object, outputAssay) <- y_normed
    
    ## return
    object
}
