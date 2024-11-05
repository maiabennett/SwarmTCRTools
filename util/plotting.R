
#########
# Plots #
#########
plotRefWeights <- function(all.weights, set.weights, directory) {
  all.long <- all.weights %>% 
    pivot_longer(cols = starts_with("CDR"), names_to = "CDR", values_to = "Weight")
  set.long <- set.weights %>% 
    pivot_longer(cols = starts_with("CDR"), names_to = "CDR", values_to = "Weight")
  epitopes <- unique(all.long$Epitope)
  ggplot(all.long, aes(x = CDR, y = Weight, color = Epitope)) +
    geom_boxplot() +
    geom_point(data =set.long, aes(x = CDR, y = Weight), color = "black") +
    facet_wrap(~Epitope) +
    scale_fill_viridis(discrete = TRUE, option = "F", begin = .3, end = .7) +
    theme(legend.position = "none") +
    labs(title = "Reference Weights",
         x = "CDR",
         y = "Weight")
}

plotScoreThresholds <- function(pr.dfs, scores.dfs, epitopes, thresholds, directory) {
  pr.long <- data.frame()
  scores.long <- data.frame()
  for (epitope in epitopes) {
    pr <- pr.dfs[[epitope]]
    pr <- pr %>% 
      mutate(Epitope = epitope)
    scores <- scores.dfs[[epitope]]
    scores <- scores %>% 
      mutate(Epitope = epitope)
    pr.long <- rbind(pr.long, pr)
    scores.long <- rbind(scores.long, scores)
  }

  (plot <- plotPrecisionRecall(pr.long, thresholds))
  print(plot)
  ggsave(paste0(directory, "thresholdPRs.png"), plot, width = 20, height = 10)

  (plot <- plotScores(scores.long, thresholds))
  print(plot)
  ggsave(paste0(directory, "scoreDistributions.png"), plot, width = 20, height = 10)

}

plotPrecisionRecall <- function(pr, thresholds) {
  labels <- c()
  for (epitope in unique(thresholds$Epitope)) {
    labels <- c(labels, paste0(epitope, 
      "\nF1: ", thresholds %>% filter(Epitope == epitope) %>% pull(F1) %>% round(3),
      #"\nPrecision: ", thresholds %>% filter(Epitope == epitope) %>% pull(Precision) %>% round(3),
      #"\nRecall: ", thresholds %>% filter(Epitope == epitope) %>% pull(Recall) %>% round(3),
      "\nThreshold: ", thresholds %>% filter(Epitope == epitope) %>% pull(Threshold) %>% round(3)))
  }
  pr <- pr %>% 
    mutate(Epitope = factor(Epitope, levels = unique(pr$Epitope), labels = labels))
  thresholds <- thresholds %>% 
    mutate(Epitope = factor(Epitope, levels = unique(thresholds$Epitope), labels = labels))
  ggplot(pr, aes(x = Recall, y = Precision, color = as.factor(Run))) +
    geom_line() +
    facet_wrap(~Epitope) +
    coord_fixed(ratio = 1) +
    scale_color_viridis(discrete = TRUE, option = "F", begin = .3, end = .7) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_point(data = thresholds, aes(x = Recall, y = Precision), color = "black", size = 5) +
    theme(legend.position = "none") 
}

plotScores <- function(scores, thresholds, fill = "Flag", labels = TRUE) {
  thresholds <- thresholds %>%
      filter(Epitope %in% unique(scores$Epitope))
  scores <- scores %>%
    na.omit()
  if (labels) {
    labels <- c()
    order <- unique(scores$Epitope)
    for (epitope in order) {
      labels <- c(labels, paste0(epitope, 
        "\nF1: ", thresholds %>% filter(Epitope == epitope) %>% pull(F1) %>% round(3),
        #"\nPrecision: ", thresholds %>% filter(Epitope == epitope) %>% pull(Precision) %>% round(3),
        #"\nRecall: ", thresholds %>% filter(Epitope == epitope) %>% pull(Recall) %>% round(3),
        "\nThreshold: ", thresholds %>% filter(Epitope == epitope) %>% pull(Threshold) %>% round(3)))
    } 
    scores <- scores %>% 
      mutate(Epitope = factor(Epitope, levels = order, labels = labels))
    thresholds <- thresholds %>%
      mutate(Epitope = factor(Epitope, levels = order, labels = labels))
  }
  ggplot(scores, aes(x = Score, fill = as.factor(!!sym(fill)))) +
    geom_density(alpha = 0.5) +
    geom_vline(data = thresholds, aes(xintercept = Threshold, linetype = "dashed"), color = "black") +
    facet_wrap(~Epitope) +
    theme(legend.position = "none") +
    scale_fill_viridis(discrete = TRUE, option = "F", begin = .3, end = .7) +
    labs(title = "Score Distributions",
         x = "Score",
         y = "Density")
}
