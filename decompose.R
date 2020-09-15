require("CAGEfightR")
require("assertthat")
require("caTools")

source("CAGEfightR_extensions/utils.R")

## Decomposes tag clusters according to pooled values and a decomposition function
## (implemented functions: summit_decompose and local_maxima_decompose)

## object: GRanges
## pooled: RangedSummarizedExperiment
## fn: function that takes as input a Views object (and optional arguments) and returns an IRanges object
## ...: arguments passed on to decomposition function
decompose <- function(object, pooled, fn=summit_decompose, ...) {

    pooled <- methods::as(rowRanges(pooled),"GRanges")

    assert_that(identical(seqlengths(object), seqlengths(pooled)))

    ## Split by strand
    message("Splitting by strand...")
    TCsByStrand <- splitByStrand(object)
    covByStrand <- splitPooled(pooled)

    ## Convert to IRangesList
    irl_plus <- methods::as(TCsByStrand$`+`,'IRangesList')
    irl_minus <- methods::as(TCsByStrand$`-`,'IRangesList')

    ## Views
    views_plus <- Views(covByStrand$`+`, irl_plus)
    views_minus <- Views(covByStrand$`-`, irl_minus)

    message("Decomposing tag clusters")
    irl_plus <- methods::as(lapply(views_plus, function(views) {cat(".");fn(views, ...)}),"IRangesList")
    cat("\n")
    irl_minus <- methods::as(lapply(views_minus, function(views) {cat(".");fn(views, ...)}),"IRangesList")
    cat("\n")

    message("Quantifying decomposed tag clusters...")
    decomposedTCs <- TCstats(coverage_plus = covByStrand$`+`,
                             coverage_minus = covByStrand$`-`,
                             tcs_plus = methods::as(irl_plus,"CompressedIRangesList"),
                             tcs_minus = methods::as(irl_minus,"CompressedIRangesList"))

    ## Carry over seqinfo and sort
    message("Preparing output...")
    seqinfo(decomposedTCs) <- seqinfo(object)
    decomposedTCs <- sort(decomposedTCs)

    ## Print some basic stats
    message("Tag clustering summary:")
    summarizeWidths(decomposedTCs)

    ## Return
    decomposedTCs
}

## Decompose tag cluster into subclusters according to CTSS expression fraction of summit CTSS expression
## Subclusters within mergeDist bp will be merged
summit_decompose <- function(views, fraction = 0.1, mergeDist=20) {

    if (length(views)==0)
        return(IRanges())

    pos <- viewApply(views, function(rle) {

        ## most common case: 1bp TC
        if (length(rle) == 1)
            return(c(1,1))

        r <- as.vector(rle)
        m <- max(r)

        k <- which(r >= fraction*m)

        ## special case: summit position only
        if (length(k) == 1)
            return(c(k,k))

        d <- diff(k)
        s <- which(d>mergeDist)

        ## no splitting of tag cluster
        if (length(s)==0)
            return(c(k[1],k[length(k)]))

        ## splitting of tag cluster
        starts <- c(k[1],k[s+1])
        ends <- c(k[s],k[length(k)])

        return(as.vector(matrix(c(starts,ends),ncol=length(starts),byrow=TRUE)))

    })

    pos <- matrix(unlist(lapply(1:length(views), function(i) {x <- pos[[i]]; lapply(seq(1,length(x),by=2), function(j) c(i,x[j],x[j+1]))})),ncol=3,byrow=TRUE)

    s <- start(views)[pos[,1]]

    ## return IRanges object
    IRanges(start=pos[,2]+s-1,end=pos[,3]+s-1)
}

## Decompose tag cluster into subclusters according to CTSS expression fraction of local maxima CTSS expression
## Performs local summit decomposition for each local maxima separately in decreasing order of expression level
## For each local summit decomposition, subclusters will be merged if within maxGap distance.
## Final subclusters within mergeDist bp will be merged
local_maxima_decompose <- function(views, fraction = 0.1, maximaDist=20, maxGap=maximaDist, mergeDist=-1) {

    if (length(views)==0)
        return(IRanges())

    pos <- viewApply(views, function(rle) {

        ## most common case: 1bp TC
        if (length(rle) == 1)
            return(c(1,1))

        r <- as.vector(rle)

        ## Local maxima
        m <- which(apply(cbind(r,runmax(r,maximaDist)),1,function(x) x[1]==x[2]))
        m <- m[order(r[m], decreasing=TRUE)]

        ## Iterate over local maxima
        starts <- ends <- c()
        for (i in m) {
            if (r[i] == 0)
                next

            k <- which(r >= fraction*r[i] & r<=r[i])

            ## special case: summit position only
            s <- i
            e <- i

            ## multiple bps
            if (length(k) > 1) {

                d <- diff(k)
                gaps <- which(d>maxGap)

                s <- c(k[1],k[gaps+1])
                e <- c(k[gaps],k[length(k)])
                idx <- max(which(s<=i))
                s <- s[idx]
                e <- e[idx]

                ## Check that sub cluster does not contain or overlap previous local maxima sub clusters
                if (any(starts %in% s:e)) {
                    es <- c(e,intersect(starts,s:e)-1)
                    e <- min(es[es>=i])
                }
                if (any(ends %in% s:e)) {
                    ss <- c(s,intersect(ends,s:e)+1)
                    s <- max(ss[ss<=i])
                }
            }

            ## set vales in sub cluster to 0 for next local maxima
            r[s:e] <- 0

            ## merge with previous cluster if proximal
            merge <- FALSE
            if (length(starts) > 0 && mergeDist>=0) {
                m <- c(which(starts-mergeDist-1 %in% s:e),which(ends+mergeDist+1 %in% s:e))
                if (length(m)>0) {
                    merge <- TRUE
                    m <- min(m)
                    starts[m] <- min(starts[m],s)
                    ends[m] <- max(ends[m],e)
                }
            }

            if (!merge) {
                ## store sub cluster
                starts <- c(starts,s)
                ends <- c(ends,e)
            }
        }

        return(as.vector(matrix(c(starts,ends),ncol=length(starts),byrow=TRUE)))
    })

    pos <- matrix(unlist(lapply(1:length(views), function(i) {x <- pos[[i]]; lapply(seq(1,length(x),by=2), function(j) c(i,x[j],x[j+1]))})),ncol=3,byrow=TRUE)

    s <- start(views)[pos[,1]]

    ## return IRanges object
    IRanges(start=pos[,2]+s-1,end=pos[,3]+s-1)
}

