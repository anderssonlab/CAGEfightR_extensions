require("CAGEfightR")
require("assertthat")

## Subsample a SummarizedExperiment to a target sequencing depth
subsampleTarget <- function(object, inputAssay = "counts", target) {

    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                is.numeric(target), target > 0)
    
    a <- assay(object,inputAssay)
    n <- ncol(a)
    nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
    s <- Matrix::colSums(a)
    d <- unlist(lapply(1:n, function(i) {
        if (s[i]<=target)
            a[nz[[i]],i]
        else
            rbinom(length(nz[[i]]),a[nz[[i]],i], target/s[i])
    }))
    keep <- which(sapply(nz,length)>0)
    assay(object, inputAssay) <- Matrix::sparseMatrix(i=unlist(nz),
                                                      j=unlist(lapply(keep, function(i) rep(i,length(nz[[i]])))),
                                                      x=d, dimnames=list(rownames(a), colnames(a)[keep]))
    
    object
}

## Subsample a SummarizedExperiment to a proportion of its sequencing depth
## Useful for saturation analyses
subsampleProportion <- function(object, inputAssay = "counts", proportion) {
 
    assert_that(methods::is(object, "SummarizedExperiment"),
                inputAssay %in% assayNames(object),
                is.numeric(proportion), proportion > 0 && proportion <=1)
    
    a <- assay(object,inputAssay)
    n <- ncol(a)
    nz <- lapply(1:n, function(i) nonzero(a[,i,drop=FALSE])[,1])
    d <- unlist(lapply(1:n, function(i) {
        rbinom(length(nz[[i]]),a[nz[[i]],i], proportion)
    }))
    keep <- which(sapply(nz,length)>0)
    assay(object, inputAssay) <- Matrix::sparseMatrix(i=unlist(nz),
                                                      j=unlist(lapply(keep, function(i) rep(i,length(nz[[i]])))),
                                                      x=d, dimnames=list(rownames(a), colnames(a)[keep]))
    
    object
}


## nonzero function from the DAPAR package (https://rdrr.io/bioc/DAPAR/src/R/utils.R)
nonzero <- function(x){
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
                      dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}
