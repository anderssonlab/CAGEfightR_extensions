require("CAGEfightR")
require("parallel")

## Object: GRanges
## ctss: SummarizedExperiment
divergentLoci <- function(object, ctss, max_gap=400, win_size=200, inputAssay="counts") {

    ctss <- methods::as(rowRanges(ctss),"GRanges")
    assert_that(checkPooled(ctss))

    message("Removing overlapping TCs by strand...")
    ## Split on strand
    TCsByStrand <- splitByStrand(object)

    ## Find overlaps between TCs on separate strands
    olaps <- findOverlaps(TCsByStrand$'-',TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)
    m_score <-  mcols(TCsByStrand$'-')$score
    p_score <-  mcols(TCsByStrand$'+')$score
    
    m_rem <- queryHits(olaps)[which(m_score[queryHits(olaps)] <= p_score[subjectHits(olaps)])]
    p_rem <- subjectHits(olaps)[which(p_score[subjectHits(olaps)] < m_score[queryHits(olaps)])]

    ## remove overlapping TCs
    TCsByStrand$'-' <- TCsByStrand$'-'[-m_rem]
    TCsByStrand$'+' <- TCsByStrand$'+'[-p_rem]

    message("Finding divergent TC pairs...")
    ## Find divergent TC pairs
    m_pad <- flank(TCsByStrand$'-', width=max_gap, start=TRUE, both=FALSE)
    pairs <- findOverlaps(m_pad,TCsByStrand$'+',maxgap=-1,type="any",select="all",ignore.strand=TRUE)

    ## Find connected components of TC pair graphs
    edge_list <- cbind(names(TCsByStrand$'-')[queryHits(pairs)],
                     names(TCsByStrand$'+')[subjectHits(pairs)])
    g <- igraph::graph_from_edgelist(edge_list,directed=FALSE)
    con <- igraph::components(g)

    message("Merging into divergent loci...")

    ## Keep only relevant TCs
    object <- object[names(con$membership)]

    ## Split TCs by loci membership
    groups <- split(object,con$membership)

    ## Merge connected components into divergent loci

    midpoint <- function(grl) {
        s <- start(grl$'+')[1]
        e <- tail(end(grl$'-'),1)

        c(s,e,round(mean(c(s,e))))
    }

    chunks <- split(1:length(groups), ceiling(10*(1:length(groups))/length(groups)))
    div_mid <- unlist(lapply(1:length(chunks), function(i) {

        r <- mclapply(groups[chunks[[i]]], function(g) {

            tcs <- splitByStrand(g)
            m <- midpoint(tcs)
            
            ## conflict?
            if (m[1] < m[2])
            {
                ## extend TCs to loci extremes
                start(tcs$'-') <- start(tcs$'-')[1]
                end(tcs$'+') <- tail(end(tcs$'+'),1)
                
                ## find overlaps between strands and prioritise by score
                olaps <- findOverlaps(tcs$'-',tcs$'+',ignore.strand=TRUE)
                ps <- tcs$'+'$score
                ms <- tcs$'-'$score
                rmm <- ms[queryHits(olaps)] < ps[subjectHits(olaps)]
                rmp <- !rmm
                if (any(rmm))
                    tcs$'-' <- tcs$'-'[-queryHits(olaps)[rmm]]
                if (any(rmp))
                    tcs$'+' <- tcs$'+'[-subjectHits(olaps)[rmp]]
                
                m <- midpoint(tcs)
            }
            
        m[3]
        },mc.cores=40)
        cat(".")
        r
    }))
        
    div_chr <- sapply(groups, function(g) as.character(seqnames(object[names(con$membership[g])[1]])))

    covByStrand <- splitPooled(ctss)
    gr <- GRanges(seqnames=div_chr,IRanges(start=div_mid,end=div_mid))
    seqinfo(gr) <- seqinfo(ctss)
    
    message("Calculating directionality...")

    win_1 <- flank(gr,width=win_size,start=TRUE,both=FALSE)
    win_2 <- flank(gr,width=win_size,start=FALSE,both=FALSE)

    ## Quantify strand-wise in flanking windows around midpoint
    M1 <- unlist(viewSums(Views(covByStrand$`-`, win_1)))
    P2 <- unlist(viewSums(Views(covByStrand$`+`, win_2)))
    M2 <- unlist(viewSums(Views(covByStrand$`-`, win_2)))
    P1 <- unlist(viewSums(Views(covByStrand$`+`, win_1)))

    ## Calculate directionality
    pooled_directionality <- (P2-M1) / (P2+M1)

    ## Test if divergent
    divergent <- (M1>P1) & (P2>M2)

    message("Calculating coverage across samples...")

    ## Quantify strand-wise in flanking windows around midpoint
    strand(win_1) <- "-"
    strand(win_2) <- "+"
    mat_2_plus <- suppressMessages(assay(quantifyClusters(ctss, win_2, inputAssay = inputAssay),inputAssay) > 0)
    mat_1_minus <- suppressMessages(assay(quantifyClusters(ctss, win_1, inputAssay = inputAssay),inputAssay) > 0)

    ## Quntify number of bidirectional cases (both strands expressed)
    bidirectional <- rowSums(mat_1_minus & mat_2_plus)

    message("Preparing output...")

    ## Build GRanges object
    start(gr) <- start(gr)-win_size
    end(gr) <- end(gr)+win_size
    gr$score <- M1+P2
    gr$thick <- div_mid
        
    mcols(gr)[, "directionality"] <- pooled_directionality
    mcols(gr)[, "bidirectionality"] <- bidirectional

    ## Remove non-divergent cases
    gr[divergent]
}

