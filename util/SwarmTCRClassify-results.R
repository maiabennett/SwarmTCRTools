
############################################
# SwarmClassify results processing methods #
############################################

processSwarmClassify <- function(results.path, exp.name, cdrs.path = NULL) {
    exp.results <- import(paste0(results.path, exp.name, "_CDRs.txt.results"))
    if (!is.null(cdrs.path)) {
        exp.cdrs <- import(paste0(cdrs.path, exp.name, "_CDRs.txt"))
        exp.results <- merge(exp.cdrs, exp.results, by = "TCR_ID")
    }
    return(exp.results)
}

processERGO <- function(results.path, exp.name) {
    exp.files <- list.files(results.path, pattern = paste0(exp.name, ".*_ERGOII-results\\.csv$"), full.names = TRUE)
    exp.results <- data.frame()
    for (file in exp.files) {
        exp.results <- rbind(exp.results, import(file) %>% select(-T.Cell.Type, -MHC))
    }
    exp.results <- exp.results %>% 
        group_by(clone.id) %>%
        filter(Score == max(Score)) %>%
        ungroup()
    return(exp.results)
}

mergePredictions <- function(swarm.path, ergo.path, exp.name, exp.path) {
    exp.sTCR.results <- import(paste0(swarm.path, exp.name, "_clonotype_results.csv"))
    exp.ergo.results <- import(paste0(ergo.path, exp.name, "_clonotype_results.csv"))
    exp.full <- import(paste0(exp.path, exp.name, "_paired_clustered.csv"))
    exp.full <- mergeSwarmClassify(exp.sTCR.results, exp.full) %>% 
        mergeERGO(exp.ergo.results, .)
    return(exp.full)
}

mergeSwarmClassify <- function(exp.results, exp.full) {
    exp.full <- merge(exp.full, 
                      (exp.results %>% 
                          mutate(sTCR.score = SCORE, sTCR.epitope = EPITOPE, sTCR.reference = REF_TCR_ID) %>%
                          select(-SCORE, -EPITOPE, -TCR_ID, -REF_TCR_ID)),
                      by = c("CDR1a", "CDR2a", "CDR3a", "CDR1b", "CDR2b", "CDR3b"), 
                      all.x = TRUE)
    return(exp.full)
}

mergeERGO <- function(exp.results, exp.full) {
    exp.full <- merge(exp.full, 
                      exp.results %>% 
                          mutate(ergo.score = Score, ergo.epitope = Peptide) %>%
                          select(-Score, -Peptide, -clone.id),                    
                      by.x = c("AV", "AJ", "BV", "BJ", "CDR3a", "CDR3b"), 
                      by.y = c("TRAV", "TRAJ", "TRBV", "TRBJ", "TRA", "TRB"), 
                      all.x = TRUE)
    return(exp.full)
}
