require("CAGEfightR")
source("CAGEfightR_extensions/utils.R")

## object: RangedSummarizedExperiment
## mask: GRanges
## mappable: GRanges

estimateNoise <- function(object, mask, mappable, map_frac=0.5, win_size=200, num_win=1e6, strand="+", inputAssay="counts", quantiles=c(0.9,0.95,0.99,0.999,0.9999,0.99999)) {

    message("Creating unmasked windows...")
    ## Create tiling windows
    genome_gr <- tileGenome(seqinfo(object),tilewidth=win_size,cut.last.tile.in.chrom=TRUE)
    ## Remove masked regions
    genome_gr <- subsetByOverlaps(genome_gr, mask, maxgap=-1, type="any", invert=TRUE)

    message("Filtering windows based on mappability...")
    hits <- findOverlaps(mappable, genome_gr, maxgap=-1, type="any")
    olaps <- pintersect(mappable[queryHits(hits)], genome_gr[subjectHits(hits)])
    df <- data.frame(subjectHits=subjectHits(hits),coverage=width(olaps) / width(genome_gr[subjectHits(hits)]))
    df <- aggregate(df, by=list(df$subjectHits), FUN=sum)

    ## Keep windows with at least map_frac fraction of bases uniquely mappable
    df <- subset(df, coverage >= map_frac)
    genome_gr <- genome_gr[df$Group.1]

    ## Sample windows randomly
    genome_gr <- genome_gr[sort(sample.int(length(genome_gr),num_win,replace=FALSE))]

    message("Quantifying expression across samples...")
    strand(genome_gr) <- strand
    expr <- suppressMessages(assay(quantifyClusters(object, genome_gr, inputAssay = inputAssay),inputAssay))

    message("Calculating statistics...")
    res <- apply(expr,2,quantile,prob=quantiles)

    res
}
