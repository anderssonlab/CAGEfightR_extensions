require("CAGEfightR")

source("CAGEfightR_extensions/utils.R")

calcComplexity <- function(object, txModels, step=1e6, CTSSunexpressed=1, geneunexpressed=9, minCTSSsupport=2) {

    object <- suppressWarnings(calcTotalTags(object, inputAssay="counts"))
    object <- suppressMessages(suppressWarnings(assignGeneID(object, geneModels = txModels, outputColumn = "geneID")))

    targets <- seq(step,max(object$totalTags),by=step)

    message("Subsampling counts")
    res <- c(list(data.frame(target=0,sample=colnames(object),totalTags=0,numberCTSSs=0,numberGenes=0)),
             bplapply(targets, function(t) {
                  
                 x <- subsampleTarget(object, "counts", t)
                 if (minCTSSsupport > 1)
                     x <- suppressMessages(subsetBySupport(x, unexpressed=CTSSunexpressed, minSamples=minCTSSsupport))

                 x <- suppressWarnings(calcTotalTags(x, inputAssay="counts"))
                 x <- suppressWarnings(calcNumberCTSSs(x, inputAssay="counts", unexpressed = CTSSunexpressed))
                 x <- suppressWarnings(calcNumberGenes(x, txModels, inputAssay="counts", unexpressed = geneunexpressed))
                 
                 df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberCTSSs=x$numberCTSSs,numberGenes=x$numberGenes)
                 df[object$totalTags<t,c("totalTags","numberCTSSs","numberGenes")] <- NA
                 
                 df
             }))
    
    do.call("rbind",res)
}
