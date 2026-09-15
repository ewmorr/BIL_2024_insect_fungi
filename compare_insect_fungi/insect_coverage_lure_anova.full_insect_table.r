##############################################################################
# Insect Sample Coverage vs. Lure -- Lineage C Full Insect Table
#
# Question: does lure identity affect trap sampling efficiency, measured as
# Good-Turing sample coverage (the estimated fraction of the true insect
# community captured by a trap-catch, "coverage" in insect_asymptotic_
# diversity.csv from insect_fungal_asymptotic_richness_iNEXT.full_insect_
# table.r)? This is a coverage/efficiency question, not a diversity question
# -- coverage is a property of how thoroughly a sample was drawn, and is the
# quantity the asymptotic (Chao) estimators use to correct for undersampling
# in the first place (see that script's header). Testing it directly against
# lure asks whether some lures draw a more complete catch (higher coverage)
# than others, independent of how diverse what they catch turns out to be.
#
# Design: one trap per site x lure combination, sampled repeatedly across the
# season (69 samples, 4 sites x 3 lures, ~6 dates each -- see project_
# organization.md's permutation-scheme conventions). Lure is therefore
# constant WITHIN trap and samples from the same trap are not independent
# replicates of "lure" -- a naive one-way ANOVA across all 69 rows would
# pseudoreplicate and ignore the site confound. Follows this project's
# established lure permutation-ANOVA pattern (identical in form to
# insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r's
# test_term_lure(), which tests the same design for the diversity metrics --
# not duplicated by sourcing that script, since this script only needs its
# already-written coverage column, not a recompute of the Chao estimates):
#   - Omnibus: coverage ~ site + lure + date vs. coverage ~ site + date,
#     nested-model F-test; null built by permuting lure labels among traps
#     WITHIN each site (999 perms), holding site/date fixed.
#   - Pairwise post-hoc (which lures differ, not just whether any do):
#     site/date-adjusted lure-coefficient contrasts pulled from the SAME
#     full 3-lure model/permutation null as the omnibus test (not a separate
#     2-lure refit -- with only one trap per site x lure, a 2-lure-only
#     permutation null is too coarse, see section 2's comment), BH-corrected
#     across the 3 pairs.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(dplyr)
library(ggplot2)
source("library/library.R")

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"
out_fig_dir <- "figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load insect asymptotic diversity table (coverage + design vars) ---

insect_asymp <- read.csv(file.path(out_data_dir, "insect_asymptotic_diversity.csv")) %>%
  mutate(date = as.Date(date))

cat(nrow(insect_asymp), "insect samples,", length(unique(insect_asymp$trap_id)), "traps,",
    length(unique(insect_asymp$site)), "sites,", length(unique(insect_asymp$lure)), "lures.\n")
print(table(insect_asymp$site, insect_asymp$lure))

n_perm <- 999

## ---- 2. Omnibus permutation ANOVA + pairwise post-hoc contrasts -----------
# Same nested-model F-test / within-site lure-permutation scheme as
# test_term_lure() in insect_fungal_asymptotic_richness_iNEXT.full_insect_
# table.r, applied to coverage instead of a diversity metric. Pairwise
# lure contrasts are extracted from this SAME full 3-lure model/permutation
# null (not a separate 2-lure refit): with only one trap per site x lure,
# restricting a permutation to 2 lures leaves just 2 traps per site to
# shuffle (2^4 = 16 possible configurations across the 4 sites), too coarse
# to resolve a real effect, whereas the full 3-lure model's null has 6^4 =
# 1296 configurations. Because site/date enter the model additively (no
# lure interaction), each lure-level coefficient IS the site/date-adjusted
# mean offset from the reference lure, so pairwise coefficient differences
# are exactly the adjusted pairwise contrasts -- computed once per
# permutation alongside the omnibus F, at no extra fitting cost.

trap_lure_map <- insect_asymp %>% distinct(trap_id, site, lure)
lures <- sort(unique(insect_asymp$lure))              # reference = lures[1]
lure_pairs <- combn(lures, 2, simplify = FALSE)

fit_full <- function(d) lm(coverage ~ site + lure + date, data = d)
fit_F_omnibus <- function(d) anova(lm(coverage ~ site + date, data = d), fit_full(d))[2, "F"]

lure_contrasts <- function(d) {
  cs <- coef(fit_full(d))
  lure_coef <- setNames(c(0, cs[paste0("lure", lures[-1])]), lures)  # reference level = 0
  vapply(lure_pairs, function(pr) unname(lure_coef[pr[1]] - lure_coef[pr[2]]), numeric(1))
}

obs_F <- fit_F_omnibus(insect_asymp)
obs_contrasts <- lure_contrasts(insect_asymp)

perm_stats <- replicate(n_perm, {
  perm_map <- trap_lure_map %>% group_by(site) %>% mutate(lure = sample(lure)) %>% ungroup()
  lure_lookup <- setNames(as.character(perm_map$lure), perm_map$trap_id)
  d2 <- insect_asymp; d2$lure <- factor(lure_lookup[d2$trap_id], levels = lures)
  c(F = fit_F_omnibus(d2), lure_contrasts(d2))
})
p_perm_omnibus <- (sum(perm_stats["F", ] >= obs_F) + 1) / (n_perm + 1)
full_aov <- anova(fit_full(insect_asymp))
r2_lure <- full_aov["lure", "Sum Sq"] / (full_aov["lure", "Sum Sq"] + full_aov["Residuals", "Sum Sq"])

omnibus_result <- data.frame(term = "lure", stat_type = "F", F = obs_F, r2 = r2_lure,
                              p_perm = p_perm_omnibus, n_perm = n_perm, n = nrow(insect_asymp))
cat("\n--- Omnibus permutation ANOVA: coverage ~ site + lure + date (lure term) ---\n")
print(omnibus_result, row.names = FALSE)
write.csv(omnibus_result, file.path(out_data_dir, "insect_coverage_lure_anova.csv"), row.names = FALSE)

pairwise_result <- bind_rows(lapply(seq_along(lure_pairs), function(i) {
  pr <- lure_pairs[[i]]
  obs_c <- obs_contrasts[i]
  perm_c <- perm_stats[i + 1, ]  # row 1 is "F"
  p_perm_pair <- (sum(abs(perm_c) >= abs(obs_c)) + 1) / (n_perm + 1)
  mean_a <- mean(insect_asymp$coverage[insect_asymp$lure == pr[1]])
  mean_b <- mean(insect_asymp$coverage[insect_asymp$lure == pr[2]])
  data.frame(lure_a = pr[1], lure_b = pr[2], mean_a = mean_a, mean_b = mean_b,
             adj_mean_diff = obs_c, p_perm = p_perm_pair)
})) %>%
  mutate(q_value = p.adjust(p_perm, method = "BH"))

cat("\n--- Pairwise lure contrasts (site+date-adjusted, BH-corrected across",
    nrow(pairwise_result), "pairs) ---\n")
print(pairwise_result, row.names = FALSE)
write.csv(pairwise_result, file.path(out_data_dir, "insect_coverage_lure_pairwise.csv"), row.names = FALSE)

## ---- 4. Boxplot: coverage by lure ------------------------------------------

p_label <- paste0("Permutation ANOVA (site+date controlled): F = ", round(obs_F, 2),
                   ", p = ", signif(p_perm_omnibus, 2))

coverage_plot <- ggplot(insect_asymp, aes(x = lure, y = coverage, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.6) +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Lure", y = "Sample coverage (Good-Turing)",
       title = "Insect trap sample coverage by lure",
       subtitle = p_label) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))
ggsave(file.path(out_fig_dir, "insect_coverage_by_lure.png"), coverage_plot, width = 5.5, height = 4.5, dpi = 150)
ggsave(file.path(out_fig_dir, "insect_coverage_by_lure.pdf"), coverage_plot, width = 5.5, height = 4.5)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
