

d <- read.csv("POLS4724_Data_Song.csv")
d$row_dist <- 10 * (d$row_num - 1)
d$col_dist <- 10 * (3 - d$col_num)

## Balance test - general
covariates <- c("w_pre", "peri_pre", "root_length_pre")
balance <- do.call(rbind, lapply(covariates, function(v) {
  m <- d[[v]][d$treatment == 1]
  d_group <- d[[v]][d$treatment == 0]
  pooled_sd <- sqrt((var(m, na.rm = TRUE) + var(d_group, na.rm = TRUE)) / 2)
  
  data.frame(
    covariate = v,
    mean_M = mean(m, na.rm = TRUE),
    mean_D = mean(d_group, na.rm = TRUE),
    SMD = (mean(m, na.rm = TRUE) - mean(d_group, na.rm = TRUE)) / pooled_sd,
    p_value = t.test(m, d_group)$p.value
  )
}))

print(balance, digits = 3, row.names = FALSE)

## Balance test - by block
design <- randomizr::declare_ra(
  blocks = d$block,
  block_m = rep(3L, length(unique(d$block)))
)

covariates <- c("w_pre", "peri_pre", "root_length_pre")
set.seed(4724)

for (v in covariates) {
  fit <- ri2::conduct_ri(
    formula = reformulate("treatment", response = v),
    assignment = "treatment",
    declaration = design,
    sharp_hypothesis = 0,
    data = d,
    sims = 100000,
    p = "two-tailed"
  )
  
  cat("\n", v, "\n", sep = "")
  print(summary(fit))
  print(plot(fit) + ggplot2::labs(title = paste("RI balance test:", v)))
}


# Results table
results <- data.frame(
  outcome = c("total_h", "max_h", "total_w", "new_w"),
  DIM = NA_real_,
  DIM_SE = NA_real_,
  DID = NA_real_,
  DID_SE = NA_real_,
  Regression = NA_real_,
  Regression_SE = NA_real_
)

## Block DIM
blocks <- sort(unique(d$block))

for (outcome in c("total_h", "max_h", "total_w", "new_w")) {
  dim_block <- var_block <- numeric(length(blocks))
  
  for (i in seq_along(blocks)) {
    in_block <- d$block == blocks[i]
    treated <- d[[outcome]][in_block & d$treatment == 1]
    control <- d[[outcome]][in_block & d$treatment == 0]
    
    dim_block[i] <- mean(treated) - mean(control)
    var_block[i] <- var(treated) / length(treated) +
      var(control) / length(control)
  }
  
  estimate <- mean(dim_block)
  SE <- sqrt(sum(var_block)) / length(blocks)
  
  i_row <- match(outcome, results$outcome)
  results$DIM[i_row] <- estimate
  results$DIM_SE[i_row] <- SE
}

## Block DID
for (outcome in c("total_w")) {
  change <- d[[outcome]] - d$w_pre
  did_block <- var_block <- numeric(length(blocks))
  
  for (i in seq_along(blocks)) {
    in_block <- d$block == blocks[i]
    treated <- change[in_block & d$treatment == 1]
    control <- change[in_block & d$treatment == 0]
    
    did_block[i] <- mean(treated) - mean(control)
    var_block[i] <- var(treated) / length(treated) +
      var(control) / length(control)
  }
  
  estimate <- mean(did_block)
  SE <- sqrt(sum(var_block)) / length(blocks)
  
  i_row <- match(outcome, results$outcome)
  results$DID[i_row] <- estimate
  results$DID_SE[i_row] <- SE
}

## Regression
for (outcome in c("total_h", "max_h", "total_w", "new_w")) {
  model <- lm(
    reformulate(
      c("treatment", "w_pre", "peri_pre", "root_length_pre",
        "row_dist", "col_dist", "factor(block)"),
      response = outcome
    ),
    data = d
  )
  
  result <- coef(summary(model))["treatment", c("Estimate", "Std. Error")]
  i_row <- match(outcome, results$outcome)
  results$Regression[i_row] <- result["Estimate"]
  results$Regression_SE[i_row] <- result["Std. Error"]
}


final_table <- data.frame(
  Outcome = results$outcome,
  `Block DIM` = sprintf("%.4f (%.4f)", results$DIM, results$DIM_SE),
  `Block DID` = sprintf("%.4f (%.4f)", results$DID, results$DID_SE),
  `Adjusted regression` = sprintf("%.4f (%.4f)", results$Regression,
                                  results$Regression_SE),
  check.names = FALSE
)

print(final_table, row.names = FALSE, right = FALSE)

# Visualization
library(ggplot2)

outcomes <- c("total_h", "max_h", "total_w", "new_w")
blocks <- sort(unique(d$block))

plot_data <- do.call(rbind, lapply(outcomes, function(outcome) {
  estimate <- se <- numeric(length(blocks))
  
  for (i in seq_along(blocks)) {
    in_block <- d$block == blocks[i]
    treated <- d[[outcome]][in_block & d$treatment == 1]
    control <- d[[outcome]][in_block & d$treatment == 0]
    
    estimate[i] <- mean(treated) - mean(control)
    se[i] <- sqrt(var(treated) / length(treated) +
                    var(control) / length(control))
  }
  
  data.frame(
    outcome = outcome,
    group = c(paste("Block", blocks), "Overall ATE"),
    estimate = c(estimate, mean(estimate)),
    se = c(se, sqrt(sum(se^2)) / length(blocks)),
    type = c(rep("Block", length(blocks)), "Overall")
  )
}))

plot_data$outcome <- factor(plot_data$outcome, levels = outcomes)
plot_data$group <- factor(
  plot_data$group,
  levels = c(rev(paste("Block", blocks)), "Overall ATE")
)
plot_data$lower <- plot_data$estimate - 1.96 * plot_data$se
plot_data$upper <- plot_data$estimate + 1.96 * plot_data$se

ggplot(plot_data, aes(estimate, group, color = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(aes(x = lower, xend = upper, yend = group),
               linewidth = 0.8, alpha = 1) +
  geom_point(aes(size = type), alpha = 1) +
  facet_wrap(~ outcome, scales = "free_x") +
  scale_color_manual(values = c(Block = "#6596C2", Overall = "#123C69")) +
  scale_size_manual(values = c(Block = 2.4, Overall = 3.5)) +
  labs(x = "Treatment effect (M minus D)", y = NULL,
       caption = "Bars show approximate 95% confidence intervals (estimate ± 1.96 SE).") +
  theme_minimal() +
  theme(legend.position = "none")

## RI
design <- randomizr::declare_ra(
  blocks = d$block,
  block_m = rep(3L, length(unique(d$block)))
)

set.seed(12345)

for (outcome in c("total_h", "max_h", "total_w", "new_w")) {
  fit <- ri2::conduct_ri(
    formula = reformulate("treatment", response = outcome),
    assignment = "treatment",
    declaration = design,
    sharp_hypothesis = 0,
    data = d,
    sims = 100000,
    p = "upper"
  )
  
  cat("\n", outcome, "\n", sep = "")
  print(summary(fit))
  print(plot(fit) + ggplot2::labs(title = paste("RI test for block DIM:", outcome)))
}







