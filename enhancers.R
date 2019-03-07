require("CAGEfightR")

## Object: GRanges
## ctss: SummarizedExperiment
divergentLoci <- function(object, ctss, max_gap=400, win_size=200, inputAssay="counts") {

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
    groups <- split(1:length(con$membership),con$membership)
    ## Merge connected components into divergent loci

    div_loci <- bplapply(groups, function(g) {
        tcs <- names(con$membership[g])
        m_tcs <- TCsByStrand$'-'[tcs[grep(";-",tcs,fixed=TRUE)]]
        p_tcs <- TCsByStrand$'+'[tcs[grep(";+",tcs,fixed=TRUE)]]

        plen <- length(p_tcs)
        mlen <- length(m_tcs)
        pi <- 1
        mi <- mlen
        
        #ms <- start(m_tcs[1])
        me <- end(m_tcs)[mi]
        ps <- start(p_tcs)[pi]
        #pe <- end(p_tcs[plen])

        if (ps < me) {
        
            m_scores <- mcols(m_tcs)$score
            p_scores <- mcols(p_tcs)$score
            
            ## Resolve conflicts
            while (ps < me) {
                mv <- m_scores[mi]
                pv <- p_scores[pi]
                
                rm_p <- mv > pv || mlen == 1
                rm_m <- pv >= mv || plen == 1
                
                if (rm_p && rm_m) {
                    if (mlen==1) rm_m <- FALSE
                    if (plen==1) rm_p <- FALSE
                }
                
                ## Remove plus strand TC and update plus strand coordinates
                if (rm_p) {
                    pi <- pi+1
                    plen <- plen-1
                    ps <- start(p_tcs)[pi]
                }
                
                ## Remove plus strand TC and update plus strand coordinates
                if (rm_m) {
                    mi <- mi-1
                    mlen <- mlen-1
                    me <- end(m_tcs)[mi]
                }
                
                ## Should never happen
                if (!rm_m && !rm_p) break
            }
        }

        list(chr=as.character(seqnames(m_tcs)[1]),mid=round(mean(c(me,ps))))
    })
    
    div_chr <- sapply(div_loci, function(x) x$chr)
    div_mid <- sapply(div_loci, function(x) x$mid)

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

