
############################################
# SwarmClassify results processing methods #
############################################

processSwarmClassify <- function(file, cdrs.path = NULL) {
    results <- read.table(file, header = TRUE)
    source <- sub(".*\\/(.*)_CDRs.*", "\\1", file)
    if (!is.null(cdrs.path)) {
        cdrs <- import(paste0(cdrs.path, source, "_CDRs.txt"))
        results <- merge(cdrs, results, by = "TCR_ID")
    }
    results <- results %>% 
        mutate(Source = source)
    return(results)
}


mergeSwarmClassify <- function(exp.results, exp.full, prefix = NULL) {
    exp.results <- exp.results %>% 
        mutate(sTCR.score = SCORE, sTCR.epitope = EPITOPE, sTCR.reference = REF_TCR_ID) %>%
        select(-SCORE, -EPITOPE, -TCR_ID, -REF_TCR_ID) 
    
    if (!is.null(prefix)) {
        exp.results <- exp.results %>% 
            rename_with(~paste0(prefix, .), starts_with("sTCR"))
    }

    exp.full <- merge(exp.full, 
        exp.results,
        by = c("CDR1a", "CDR2a", "CDR3a", "CDR1b", "CDR2b", "CDR3b"), 
        all.x = TRUE)
    return(exp.full)
}

mergeTCRDist <- function(exp.results, exp.full, prefix = NULL) {
    exp.results <- exp.results %>% 
        mutate(TCRdist.score = SCORE, TCRdist.epitope = EPITOPE, TCRdist.reference = REF_TCR_ID) %>%
        select(-SCORE, -EPITOPE, -TCR_ID, -REF_TCR_ID) 

    if (!is.null(prefix)) {
        exp.results <- exp.results %>% 
            rename_with(~paste0(prefix, .), starts_with("TCRdist"))
    }

    exp.full <- merge(exp.full, 
        exp.results,
        by = c("CDR1a", "CDR2a", "CDR3a", "CDR1b", "CDR2b", "CDR3b"), 
        all.x = TRUE)
    return(exp.full)
}


thresholdSwarmClassify <- function(exp.results, thresholds) {
    epitopes <- unique(thresholds$Epitope)
    results <- data.frame()
    for (epitope in epitopes) {
        threshold <- thresholds %>% 
            filter(Epitope == epitope) %>% 
            select(Threshold) %>% 
            pull()
        epitope.results <- exp.results %>% 
            filter(EPITOPE == epitope, SCORE >= threshold)
        results <- rbind(results, epitope.results)
    }
    return(exp.results)
}

thresholdByDist <- function(exp.results, threshold) {
    exp.results <- exp.results %>% 
        filter(SCORE >= threshold)
    return(exp.results)
}
