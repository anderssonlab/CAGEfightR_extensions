require("CAGEfightR")
require("assertthat")
require("caTools")

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

local_maxima_decompose <- function(views, fraction = 0.1, maximaDist=20) {

    if (length(views)==0)
        return(IRanges())

    pos <- viewApply(views, function(rle) {
        
        ## most common case: 1bp TC
        if (length(rle) == 1)
            return(c(1,1))
        
        r <- as.vector(rle)

        ## Local maxima
        m <- which(apply(cbind(r,runmax(r,maximaDist*2+1)),1,function(x) x[1]==x[2]))
        m <- m[order(r[m], decreasing=TRUE)]

        ## Iterate over local maxima
        starts <- ends <- c()
        maxima <- c()
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
                gaps <- which(d>maximaDist/2)
                
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

            ## store sub cluster
            starts <- c(starts,s)
            ends <- c(ends,e)
            maxima <- c(maxima,i)
        }

        return(as.vector(matrix(c(starts,ends),ncol=length(starts),byrow=TRUE)))
    })

    pos <- matrix(unlist(lapply(1:length(views), function(i) {x <- pos[[i]]; lapply(seq(1,length(x),by=2), function(j) c(i,x[j],x[j+1]))})),ncol=3,byrow=TRUE)

    s <- start(views)[pos[,1]]

    ## return IRanges object
    IRanges(start=pos[,2]+s-1,end=pos[,3]+s-1)
}

### Helper functions not exported by CAGEfightR

TCstats <- function(coverage_plus, coverage_minus, tcs_plus, tcs_minus) {
                                        # Check classes
    stopifnot(methods::is(coverage_plus, "SimpleRleList"),
              methods::is(coverage_minus, "SimpleRleList"),
              methods::is(tcs_plus, "CompressedIRangesList"),
              methods::is(tcs_minus, "CompressedIRangesList"))

                                        # Check seqlevels
    stopifnot(length(unique(list(names(coverage_plus),
                                 names(tcs_plus),
                                 names(coverage_minus),
                                 names(tcs_minus)))) == 1)

                                        # Obtain views
    views_plus <- Views(coverage_plus, tcs_plus)
    views_minus <- Views(coverage_minus, tcs_minus)

                                        # Calculate Sums
    sum_plus <- unlist(viewSums(views_plus))
    sum_minus <- unlist(viewSums(views_minus))

                                        # Find peaks
    ranges_plus <- viewRangeMaxs(views_plus)
    ranges_minus <- viewRangeMaxs(views_minus)
    ranges_plus <- resize(unlist(ranges_plus), width = 1, fix = "center")
    ranges_minus <- resize(unlist(ranges_minus), width = 1, fix = "center")

                                        # Merge into GRanges
    TCs <- c(GRanges(tcs_plus,
                     strand = "+",
                     score = sum_plus,
                     thick = ranges_plus),
             GRanges(tcs_minus,
                     strand = "-",
                     score = sum_minus,
                     thick = ranges_minus))

                                        # Names as IDs for both ranges and peaks
    TC_ids <- paste0(seqnames(TCs), ":", start(TCs), "-", end(TCs), ";", strand(TCs))
    names(TCs) <- TC_ids
    names(TCs$thick) <- TC_ids

                                        # Return
    TCs
}

summarizeWidths <- function(gr) {
                                        # Checks
    stopifnot(methods::is(gr, "GRanges"))

                                        # Cut up widths
    x <- cut(width(gr), breaks = c(1, 10, 100, 1000, Inf), labels = c(">=1", ">=10",
                                                                      ">=100", ">=1000"), include.lowest = TRUE)

                                        # Get freqs and props
    y <- table(Width = x)
    z <- prop.table(y)

                                        # Format to data.frame
    w <- merge(as.data.frame(y, responseName = "Count"),
               as.data.frame(z, responseName = "Percent"))

                                        # Add Total row
    w <- rbind(data.frame(Width = "Total",
                          Count = sum(w$Count),
                          Percent = sum(w$Percent)), w)

                                        # Reformat to percent
    w$Percent <- paste0(format(w$Percent * 100, digits = 1), " %")

                                        # To string and message
    s <- paste(utils::capture.output(print(w, row.names = FALSE)), collapse = "\n")
    message(s)
    }

splitByStrand <- function(object) {
    split(object, strand(object))
}

splitPooled <- function(object){

    ## Split by strand
    o <- splitByStrand(object)

    ## Calculate coverage
    o <- lapply(o, coverage, weight="score")

    ## Round to handle floating point errors
    o <- lapply(o, round, digits=9)

    ## Return
    o
}
