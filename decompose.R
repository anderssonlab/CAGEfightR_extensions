require("CAGEfightR")
require("assertthat")
require("caTools")
require("bcp")

source("CAGEfightR_extensions/utils.R")

## Decompose tag clusters according to sample correlations
## object: GRanges
## ctss: RangedSummarizedExperiment
decomposeCorr <- function(object, ctss, fn=corr_decompose, ...) {

    assert_that(identical(seqlengths(object), seqlengths(ctss)))

    ## Split by strand
    message("Splitting and intersecting data...")
    TCsByStrand <- splitByStrand(object)
    covByStrand <- splitPooled(methods::as(rowRanges(ctss),"GRanges"))
    ctssByStrand <- splitByStrand(ctss)

    ## Convert TCs to GRangesList split by chromosome
    grl_plus <- methods::as(split(TCsByStrand$`+`, seqnames(TCsByStrand$`+`)),'GRangesList')
    grl_minus <- methods::as(split(TCsByStrand$`-`, seqnames(TCsByStrand$`-`)),'GRangesList')

    ## Create CTSS intersects
    ctss_plus <- lapply(grl_plus, function(gr) subsetByOverlaps(ctssByStrand$`+`, gr))
    ctss_minus <- lapply(grl_minus, function(gr) subsetByOverlaps(ctssByStrand$`-`, gr))

    message("Decomposing tag clusters")
    irl_plus <- methods::as(lapply(seq_along(ctss_plus), function(i) fn(ctss_plus[[i]], grl_plus[[i]], ...)),"IRangesList")
    irl_minus <- methods::as(lapply(seq_along(ctss_minus), function(i) fn(ctss_minus[[i]], grl_minus[[i]], ...)),"IRangesList")

    names(irl_plus) <- names(grl_plus)
    names(irl_minus) <- names(grl_minus)

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
    irl_plus <- methods::as(lapply(views_plus, function(views) fn(views, ...)),"IRangesList")
    irl_minus <- methods::as(lapply(views_minus, function(views) fn(views, ...)),"IRangesList")

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


corr_decompose <- function(rse, gr, assay="TPM", thres=0.25, merge=TRUE, scale=TRUE) {

    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

    if (length(rse)==0)
        return(IRanges())

    fo <- findOverlaps(rse, gr)
    hits <- sort(unique(subjectHits(fo)))

    assert_that(length(hits) == length(gr))

    pos <- bplapply(hits, function(i) {

        q <- queryHits(fo)[which(subjectHits(fo) == i)]

        ## most common case: 1bp TC
        if (length(q) == 1)
            return(c(1,1))

        g <- gr[i]
        r <- rse[q]

        mat <- as.matrix(t(assay(r,assay)))
        if (scale)
            mat <- as.matrix(apply(mat,2,scale))
        corr <- cor(mat,method="pearson")
        eig <- eigen(corr)$vectors[,1]
        bcp_eig <- bcp(eig)

        bp <- which(bcp_eig$posterior.prob > thres)
        rel.pos <- start(rowRanges(r)) - start(rowRanges(r))[1] + 1

        if (length(bp) == 0)
            return(c(1,rel.pos[length(q)]))

        if (max(bp) == length(q))
            bp <- bp[-length(bp)]

        if (length(bp) == 0)
            return(c(1,rel.pos[length(q)]))

        spl <- splitAt(1:length(q), bp+1)

        starts <- sapply(spl, function(x) x[1])
        ends <- sapply(spl, function(x) x[length(x)])

        if (merge) {
            nb.corr <- sapply(1:(length(spl)-1), function(j) mean(corr[spl[[j]],spl[[j+1]]]))
            keep <- 1:(length(starts)-1)

            ## Merge decomposed clusters if positively correlated
            if (any(nb.corr > 0))
                keep <- which(nb.corr < 0)

            starts <- c(starts[1],starts[keep+1])
            ends <- c(ends[keep],ends[length(ends)])
        }

        return(as.vector(matrix(c(rel.pos[starts],rel.pos[ends]),
                                ncol=length(starts),byrow=TRUE)))
    })

    pos <- matrix(unlist(lapply(hits, function(i) {
        x <- pos[[i]]
        lapply(seq(1,length(x),by=2), function(j) c(i,x[j],x[j+1]))
    })),ncol=3,byrow=TRUE)

    s <- start(gr)[pos[,1]]

    ## return IRanges object
    IRanges(start=pos[,2]+s-1,end=pos[,3]+s-1)
}

## Decompose tag cluster into subclusters according to CTSS expression fraction of local maxima CTSS expression
## Performs local summit decomposition for each local maxima separately in decreasing order of expression level
## For each local summit decomposition, subclusters will be merged if within maxGap distance.
## If smoothPad>0, neighbouring non-zero CTSSs within smoothPad distance of CTSSs fulfilling the summit fraction criterion will also be included.
## Final subclusters within mergeDist bp will be merged
local_maxima_decompose <- function(views, fraction = 0.1, maximaDist=20, maxGap=maximaDist, mergeDist=-1, smoothPad=0) {

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

            if (smoothPad>0) {
                nk <- k
                for (j in 1:smoothPad) {
                    nk <- c(nk,k+j,k-j)
                }
                k <- sort(unique(nk))
                if (any(k<1))
                    k <- k[-which(k<1)]
                if (any(k>length(r)))
                    k <- k[-which(k>length(r))]
                if (any(r[k]==0))
                    k <- k[-which(r[k]==0)]
            }

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

            starts <- c(starts,s)
            ends <- c(ends,e)
        }

        o <- order(starts)
        starts <- starts[o]
        ends <- ends[o]

        ## merge clusters if proximal
        if (length(starts) > 0 && mergeDist>=0) {
            d <- starts[-1]-ends[-length(ends)]
            gaps <- which(d>mergeDist)
            starts <- c(starts[1],starts[gaps+1])
            ends <- c(ends[gaps],ends[length(ends)])
        }

        return(as.vector(matrix(c(starts,ends),ncol=length(starts),byrow=TRUE)))
    }, simplify=FALSE)

    pos <- matrix(unlist(lapply(1:length(views), function(i) {x <- pos[[i]]; lapply(seq(1,length(x),by=2), function(j) c(i,x[j],x[j+1]))})),ncol=3,byrow=TRUE)

    s <- start(views)[pos[,1]]

    ## return IRanges object
    IRanges(start=pos[,2]+s-1,end=pos[,3]+s-1)
}

