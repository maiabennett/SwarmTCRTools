
##########################################
# SwarmTCR training + validation methods #
##########################################

## Reference set preparation ##

# There are three files per epitope required to train and validate SwarmTCR: 
# a reference file (distinct from train and test) which contains only positive binding pairs (flag 1),
# and training and test files which contains positive and negative binding pairs (flag 1 and 0, respectively)
splitReferenceData <- function(pos.data, neg.data = NULL, target.epitopes, directory = "./reference-data/", ref.prop = 0.3, train.prop = 0.7, test.prop = 0.3, n.folds = 25, neg.prop = NULL) {

    pos.data <- pos.data %>% 
        formatSwarmTCR(keep_epitopes = TRUE)


    # Create a list of epitope-specific dataframes for each epitope in epitopes
    epitope.dfs <- list()
    for (epitope in target.epitopes) {
        epitope.dfs[[epitope]] <- pos.data %>% 
        filterEpitope(epitope) %>%
        select(-Epitope)
    }

    # Split the epitope reference datasets into train, test, and ref sets
    for (epitope in target.epitopes) {
        epitope.df <- epitope.dfs[[epitope]] 

        for (i in 1:n.folds) {
        # Assign a portion of the epitope-specific data to the train/test reference set
        ref.data <- epitope.df %>% 
            dplyr::slice_sample(n = floor(ref.prop * nrow(epitope.df))) %>% 
            mutate(Flag = 1) %>%
            select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)

        # If a separate negative dataset is not provided, simply split the positive dataset
        if (is.null(neg.data)) {
            # Assign the remaining portion of the epitope-specific data to the train/test reference set with flag 1
            # Assign all remaining non-specific reference data to the train/test reference set with flag 0
            train.test.data <- pos.data %>%
                filter(!str_detect(Epitope, epitope)) %>%
                mutate(Flag = 0) %>%
                select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)

            # Optionally, limit the number of negative pairs to a proportion of the positive pairs
            if (!is.null(neg.prop)) {
                num.pos <- nrow(ref.data) * (1 - ref.prop)
                num.neg <- floor(neg.prop * num.pos)
                
                # Check to ensure the target number of negative pairs is less than the number of negative pairs available, then randomly sample the negative pairs to target proportion
                if (num.neg < nrow(train.test.data)) {
                    train.test.data <- train.test.data %>% 
                        dplyr::slice_sample(n = num.neg)
                }
            }

            train.test.data <- train.test.data %>% 
                bind_rows(pos.data %>%
                    filter(str_detect(Epitope, epitope)) %>%
                    mutate(Flag = 1) %>%
                    anti_join(ref.data, by = "TCR_ID") %>%
                    select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
                    )

            # Split the train/test reference set into train and test sets, keeping the same proportion of flag 1 and 0 pairs
            train.data <- train.test.data %>% 
                dplyr::slice(createDataPartition(train.test.data$Flag, p = train.prop, list = FALSE))
            test.data <- train.test.data %>%
                anti_join(train.data, by = "TCR_ID")
        } else {
            # If negative data is provided, train/test data will be made up of epitope-specific references from the positive data and non epitope-specific references from the negative data
            neg.data <- neg.data %>% 
                formatSwarmTCR(keep_epitopes = TRUE)
            # Assign the remaining portion of the epitope-specific data to the train/test reference set with flag 1
            # Assign all remaining non-specific reference data to the train/test reference set with flag 0
            train.test.data <- neg.data %>%
                filter(!str_detect(Epitope, epitope)) %>%
                mutate(Flag = 0) %>%
                select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
            
            # Optionally, limit the number of negative pairs to a proportion of the positive pairs
            if (!is.null(neg.prop)) {
                num.pos <- nrow(ref.data) * (1 - ref.prop)
                num.neg <- floor(neg.prop * num.pos)
                
                # Check to ensure the target number of negative pairs is less than the number of negative pairs available, then randomly sample the negative pairs to target proportion
                if (num.neg < nrow(train.test.data)) {
                    train.test.data <- train.test.data %>% 
                        dplyr::slice_sample(n = num.neg)
                }
            }

            train.test.data <- train.test.data %>% 
                bind_rows(pos.data %>%
                    filter(str_detect(Epitope, epitope)) %>%
                    mutate(Flag = 1) %>%
                    anti_join(ref.data, by = "TCR_ID") %>%
                    select(TCR_ID, Flag, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b)
                    )

            # Split the train/test reference set into train and test sets, keeping the same proportion of flag 1 and 0 pairs
            train.data <- train.test.data %>% 
                dplyr::slice(createDataPartition(train.test.data$Flag, p = train.prop, list = FALSE))
            test.data <- train.test.data %>%
                anti_join(train.data, by = "TCR_ID")
        }
        
        dir.create(paste0(directory, epitope, "/fold", i, "/"), recursive = TRUE, showWarnings = FALSE)
        write_delim(ref.data, paste0(directory, epitope, "/fold", i, "/", "Reference.txt"))
        write_delim(train.data, paste0(directory, epitope, "/fold", i, "/", "Train.txt"))
        write_delim(test.data, paste0(directory, epitope, "/fold", i, "/", "Test.txt"))
        }
    }
}


# To use SwarmTCRClassify after training, each epitope needs a reference set containing all positive binding pairs (no split)
# Or, use the pre-trained weights and references provided in SwarmClassify
makeClassifierReference <- function(data, epitopes, directory = "./reference-data/", epitope.directory = TRUE) {
    for (epitope in epitopes) {
        epitope.df <- data %>% 
        filterEpitope(epitope) %>%
        select(-Epitope) %>% 
        formatSwarmTCR(flag = NULL, keep_epitopes = FALSE)
        if (epitope.directory) {
        dir.create(paste0(directory, epitope, "/"), recursive = TRUE, showWarnings = FALSE)
        write_delim(epitope.df, paste0(directory, epitope, "/", epitope, "_Classifier-Reference.txt"))
        } else {
        write_delim(epitope.df, paste0(directory, epitope, "_Classifier-Reference.txt"))
        }
    }
}