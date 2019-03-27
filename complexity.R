require("CAGEfightR")

source("CAGEfightR_extensions/utils.R")

calcComplexity <- function(object, inputAssay = "counts", step=1e6, unexpressed=0) {

    object <- calcTotalTags(object, inputAssay=inputAssay)

    targets <- seq(step,max(object$totalTags),by=step)

    res <- lapply(targets, function(t) {
        message("subsampling to ",t,"...")

        x <- subsampleTarget(object, inputAssay, t)
        x <- calcNumberCTSSs(x, inputAssay, unexpressed = unexpressed)
        x <- calcTotalTags(x, inputAssay=inputAssay)

        df <- data.frame(target=t,sample=colnames(x),totalTags=x$totalTags,numberCTSSs=x$numberCTSSs)
        df[object$totalTags<t,c("totalTags","numberCTSSs")] <- NA

        df
    })
    
    do.call("rbind",res)
}
