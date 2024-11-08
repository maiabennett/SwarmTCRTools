
#########
# Plots #
#########
plotRefWeights <- function(all.weights, set.weights) {
  all.long <- all.weights %>% 
    pivot_longer(cols = starts_with("CDR"), names_to = "CDR", values_to = "Weight")
  set.long <- set.weights %>% 
    pivot_longer(cols = starts_with("CDR"), names_to = "CDR", values_to = "Weight")
  epitopes <- unique(all.long$Epitope)
  plot <- ggplot(all.long, aes(x = CDR, y = Weight, color = Epitope)) +
    geom_boxplot() +
    geom_point(data =set.long, aes(x = CDR, y = Weight), color = "black") +
    facet_wrap(~Epitope) +
    scale_fill_viridis(discrete = TRUE, option = "F", begin = .3, end = .7) +
    theme(legend.position = "none") +
    labs(title = "Reference Weights",
         x = "CDR",
         y = "Weight")

  return(plot)
}

plotPrecisionRecall <- function(pr.dfs, epitopes, thresholds = NULL) {

  pr.long <- data.frame()
  for (epitope in epitopes) {
    pr <- pr.dfs[[epitope]]
    pr <- pr %>% 
      mutate(Epitope = epitope)
    pr.long <- rbind(pr.long, pr)
  }

  average.pr <- pr.long %>%
    group_by(Epitope, Recall) %>%
    summarize(Precision = mean(Precision)) %>%
    ungroup()

  pr.auc = data.frame()
  for (epitope in epitopes) {
    pr <- pr.long %>% filter(Epitope == epitope)
    ord <- order(pr$Recall)
    recall <- pr$Recall[ord]
    precision <- pr$Precision[ord]
    pr.auc <- rbind(pr.auc, data.frame(
        Epitope = epitope, 
        PR_AUC = sum(diff(recall) * (precision[-length(precision)] + precision[-1]) / 2)))
  }

  average.pr <- average.pr %>% 
    mutate(Epitope = factor(Epitope, levels = unique(average.pr$Epitope), labels = paste0(unique(average.pr$Epitope), "\nPR AUC: ", pr.auc$PR_AUC %>% round(3))))

  plot <- ggplot(average.pr, aes(x = Recall, y = Precision, color = Epitope)) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = "Average Precision-Recall Curve", x = "Recall", y = "Precision") +
    facet_wrap(~Epitope) +
    theme_minimal() +
    theme(legend.position = "none")

  if (!is.null(thresholds)) {
    thresholds <- thresholds %>% 
      left_join(pr.auc, by = "Epitope") %>% 
      mutate(Epitope = paste0(Epitope, "\nPR AUC: ", PR_AUC %>% round(3))) %>% 
      mutate(label = paste0("F1: ", F1 %>% round(3), "\nFbeta: ", Fbeta %>% round(3), "\nPrecision: ", Precision %>% round(3), "\nRecall: ", Recall %>% round(3), "\nThreshold: ", Threshold %>% round(3)))
    plot <- plot +
      geom_point(data = thresholds, aes(x = Recall, y = Precision), color = "black", size = 5) +
      geom_text(data = thresholds, aes(x = Recall, y = Precision, label = label), color = "black", hjust = 0, vjust = 0)
  }

  return(plot)
}

comparePrecisionRecall <- function(pr.dfs, epitopes, color = "Method") {
  pr.long <- data.frame()
  for (epitope in epitopes) {
    pr <- pr.dfs[[epitope]]
    pr <- pr %>% 
      mutate(Epitope = epitope)
    pr.long <- rbind(pr.long, pr)
  }

  average.pr <- pr.long %>%
    group_by(Epitope, Recall, !!sym(color)) %>%
    summarize(Precision = mean(Precision)) %>%
    ungroup()

  color.values <- unique(pr.long[[color]])

  pr.auc = data.frame()
  for (epitope in epitopes) {
    for (value in color.values) {
      pr <- pr.long %>% filter(Epitope == epitope, !!sym(color) == value)
      ord <- order(pr$Recall)
      recall <- pr$Recall[ord]
      precision <- pr$Precision[ord]
      pr.auc <- rbind(pr.auc, data.frame(
          Epitope = epitope, 
          Method = value,
          PR_AUC = sum(diff(recall) * (precision[-length(precision)] + precision[-1]) / 2)))
    }
  }

    generate.labels <- function(epitope, color.values, pr.auc, color) {
    labels <- sapply(color.values, function(value) {
      auc_value <- pr.auc %>% filter(Epitope == epitope, !!sym(color) == value) %>% pull(PR_AUC) %>% round(3)
      paste0(value, " PR AUC: ", auc_value)
    }, USE.NAMES = FALSE)
    return(paste(labels, collapse = "\n"))
  }

  average.pr <- average.pr %>% 
    mutate(Epitope = factor(Epitope, levels = unique(average.pr$Epitope), 
      labels = sapply(unique(average.pr$Epitope), function(epitope) {
        paste0(epitope, "\n", generate.labels(epitope, color.values, pr.auc, color))
      })))

  plot <- ggplot(average.pr, aes(x = Recall, y = Precision, color = !!sym(color))) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = "Average Precision-Recall Curve", x = "Recall", y = "Precision") +
    facet_wrap(~Epitope) +
    theme_minimal()  


  if (color == "Epitope") {
    plot <- plot +
      theme(legend.position = "none")
  }

  return(plot)
}

plotScores <- function(scores.dfs, epitopes, thresholds = NULL, fill = "Flag", labels = TRUE) {

  scores.long <- data.frame()
  for (epitope in epitopes) {
    scores <- scores.dfs[[epitope]]
    scores <- scores %>% 
      mutate(Epitope = epitope)
    scores.long <- rbind(scores.long, scores)
  }

  
  if (!is.null(thresholds)) {
    thresholds <- thresholds %>%
        filter(Epitope %in% unique(scores.long$Epitope))
    scores.long <- scores.long %>%
      na.omit()
    if (labels) {
      labels <- c()
      order <- unique(scores.long$Epitope)
      for (epitope in order) {
        labels <- c(labels, paste0(epitope, 
          "\nThreshold: ", thresholds %>% filter(Epitope == epitope) %>% pull(Threshold) %>% round(3)))
      } 
      scores.long <- scores.long %>% 
        mutate(Epitope = factor(Epitope, levels = order, labels = labels))
      thresholds <- thresholds %>%
        mutate(Epitope = factor(Epitope, levels = order, labels = labels))
    }

  plot <- ggplot(scores.long, aes(x = Score, fill = as.factor(!!sym(fill)))) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Epitope, scales = "free") +
    # theme(legend.position = "none") +
    scale_fill_viridis(discrete = TRUE, option = "F", begin = .3, end = .7) +
    labs(title = "Score Distributions",
         x = "Score",
         y = "Density")

    if (!is.null(thresholds)) {
    plot <- plot +
      geom_vline(data = thresholds, aes(xintercept = Threshold, linetype = "dashed"), color = "black")
    }
  }

  return(plot)
}

plotPRAUC <- function(avg.pr.auc, per.fold.pr.auc, fill = "Epitope") {

  plot <- ggplot(per.fold.pr.auc, aes(x = Epitope, y = PR_AUC, fill = !!sym(fill))) +
    geom_boxplot()  +
    labs(title = "PR AUC",
         x = "Epitope",
         y = "PR AUC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (fill == "Epitope") {
    plot <- plot +
      theme(legend.position = "none") +
    geom_point(data = avg.pr.auc, aes(x = Epitope, y = PR_AUC), color = "black", shape = 8)
  }

  return(plot)

}
