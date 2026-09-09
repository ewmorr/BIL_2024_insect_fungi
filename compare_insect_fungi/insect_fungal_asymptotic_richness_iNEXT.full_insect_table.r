##############################################################################
# Insect vs. Fungal Asymptotic (Extrapolated) Diversity -- iNEXT, Full Insect
# Table Variant
#
# Extends insect_fungal_alpha_diversity.full_insect_table.r's question --
# does insect diversity track fungal diversity, and do site/lure/date affect
# each community's own diversity -- to EXTRAPOLATED (asymptotic) diversity
# instead of observed diversity. Observed richness/Shannon/Simpson undercount
# true diversity whenever a community isn't exhaustively sampled; the iNEXT
# package's Chao-type estimators use the shape of each sample's rare-species
# tail (singletons/doubletons) to estimate the richness/diversity that would
# be seen with infinite additional sampling effort ("asymptotic" diversity).
#
# Insect table: all families, RAW counts -- NO singleton-taxon filter, NO
# >=5-sample prevalence filter. This departs from insect_fungal_alpha_
# diversity.full_insect_table.r and Lineage C generally, which apply a global
# singleton filter (colSums > 1) as noise-floor cleanup -- see project_
# organization.md's Lineage C section and "When to prevalence-filter vs.
# not". That cleanup is WRONG for Chao-type estimators specifically: Chao1
# richness (S_obs + f1^2/(2*f2)), Chao-Shannon and Chao-Simpson are all
# functions of the sample's rare tail -- the singleton count f1 (entering
# squared for richness) and doubleton count f2, plus the Good-Turing
# coverage C_hat = 1 - f1/n. A global singleton (colSums == 1) is by
# definition a single individual in a single sample, i.e. a within-sample
# singleton, so filtering these strips f1 mass directly and unevenly across
# samples: it inflates estimated coverage and biases Chao richness/Shannon
# downward (measured: up to ~2x on the worst-sampled insect samples).
# Asymptotic estimation IS the method that reads the rare tail as signal, so
# it must see the untouched counts -- "prevalence-filtering out rare taxa
# before asking how many rare taxa we are missing would be circular", and
# the singleton filter is the same circularity in milder form. The fungal
# side was already raw (Kingdom==Fungi only, no singleton/prevalence
# filter), so this brings the insect side into line. (Switched 2026-09-09
# after the filtered version was found to distort the estimates; see
# iterative_analysis_updates.md.)
#
# Fungal table: Kingdom==Fungi filtered, but UNRAREFIED (raw ASV counts) --
# a deliberate departure from every other alpha-diversity script in this
# project, which rarefies-and-averages (multiple_subsamples(), depth 5000 x
# 100 iterations) to control for uneven sequencing depth before computing
# observed metrics. That control is unnecessary here: asymptotic estimation
# IS a principled correction for uneven sampling effort (it's the whole
# point of extrapolating each sample's own accumulation curve to its own
# asymptote), and rarefying first would throw away exactly the deep-sample
# information (the shape of the rare tail beyond depth 5000) the estimator
# needs to extrapolate accurately, and artificially cap "observed" richness
# at whatever the rarefaction depth happens to produce.
#
# iNEXT internals note: iNEXT::iNEXT() computes full sample-size-based
# rarefaction/extrapolation curves (default 40 knots) plus a bootstrapped
# asymptotic table ($AsyEst) in one call, but the curve computation over
# ~140 samples (fungal samples up to ~800k reads) does not finish in
# reasonable time (a single largest-fungal-sample call did not return in
# 120s). iNEXT also exports the lighter-weight ChaoRichness()/ChaoShannon()/
# ChaoSimpson() -- the same underlying Chao-type point/SE/CI estimators
# WITHOUT the curve machinery (~1-2s even for the largest fungal sample here,
# verified to reproduce iNEXT()'s $AsyEst point estimates exactly on a test
# subset) -- used throughout below instead.
#
# Two known iNEXT quirks worked around here:
# - ChaoShannon()/ChaoSimpson() default to Shannon entropy / Gini-Simpson
#   index (transform = FALSE); transform = TRUE returns the Hill-number
#   ("effective number of species") form instead, which is what iNEXT()'s
#   $AsyEst reports and what's used throughout below.
# - ChaoRichness()/ChaoSimpson() name their SE column "Est_s.e." (trailing
#   period); ChaoShannon() names the same column "Est_s.e" (no period) -- an
#   inconsistency in the package itself. Extracted by COLUMN POSITION below
#   (chao_row()), never by name, to avoid silently returning NA.
#
# NOTE on naming: the "simpson" metric here is Hill number q=2 (inverse
# Simpson DIVERSITY, higher = more even/diverse) -- the OPPOSITE direction
# from "simpson_dominance" (D = sum(p_i^2), higher = more dominated by few
# taxa) used in insect_fungal_alpha_diversity*.r elsewhere in this project.
# Don't compare the two scripts' "simpson" numbers directly.
#
# Scope: the cross-community correlation (the user's primary question --
# does extrapolated insect richness correlate with extrapolated fungal
# richness) and the observed-vs-asymptotic comparison report BOTH observed
# and asymptotic values side by side. The heavier permutation-based
# site/lure/date covariate tests and the lure-interaction test -- mirroring
# insect_fungal_alpha_diversity.full_insect_table.r's sections 5b/6, the
# "same type of statistical analysis pipeline as the raw estimates" the user
# asked for -- are run on the ASYMPTOTIC values only (the new content this
# script adds); observed-value covariate significance already exists in that
# script (under different metric definitions -- raw Shannon/Simpson
# dominance, rarefied-and-averaged fungal side -- so isn't just re-derived
# here a second time).
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(iNEXT)
source("library/library.R")

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"
out_fig_dir <- "figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data (full insect table RAW counts, all families, raw fungal counts) ----
# Matching logic follows insect_fungal_alpha_diversity.full_insect_table.r's
# section 1, with two deliberate departures (both in the header): the fungal
# table is left unrarefied, and the insect table keeps its raw counts with NO
# singleton filter (Chao estimators key on the singleton/doubleton tail).

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

cat("Full insect table:", nrow(sp_tab), "taxa across", length(unique(sp_tab$Family)), "families.\n")
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

# NO colSums > 1 singleton filter here -- raw counts feed the Chao estimators (see header).
insect_full <- sp_tab.t[rowSums(sp_tab.t) > 0, ]
cat("After dropping empty samples only (no singleton filter):", nrow(insect_full), "samples,",
    ncol(insect_full), "insect taxa.\n")

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
cat(sum(is_fungus), "of", ncol(fungal), "ASVs confirmed Kingdom == k__Fungi (dropping",
    sum(!is_fungus), "non-fungal/unidentified ASVs).\n")
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect), "insect taxa (full table),",
    ncol(fungal), "fungal ASVs (unrarefied).\n")

## ---- 2. Asymptotic-diversity helper: Chao richness + Hill-number Shannon/Simpson ----

chao_row <- function(x) unname(unlist(x[1, 1:5]))  # Observed, Estimator, SE, LCL, UCL -- by position, see header

compute_asymptotic_metrics <- function(count_matrix) {
  ids <- rownames(count_matrix)
  rows <- lapply(ids, function(id) {
    x <- as.numeric(count_matrix[id, ])
    info <- DataInfo(data.frame(v = x), datatype = "abundance")
    rich <- chao_row(suppressWarnings(suppressMessages(ChaoRichness(x))))
    sh   <- chao_row(suppressWarnings(suppressMessages(ChaoShannon(x, transform = TRUE))))
    si   <- chao_row(suppressWarnings(suppressMessages(ChaoSimpson(x, transform = TRUE))))
    data.frame(
      sample_id = id, n_individuals = info$n[1], coverage = info$SC[1],
      richness_observed = rich[1], richness_asymptotic = rich[2],
      richness_se = rich[3], richness_lcl = rich[4], richness_ucl = rich[5],
      shannon_observed = sh[1], shannon_asymptotic = sh[2],
      shannon_se = sh[3], shannon_lcl = sh[4], shannon_ucl = sh[5],
      simpson_observed = si[1], simpson_asymptotic = si[2],
      simpson_se = si[3], simpson_lcl = si[4], simpson_ucl = si[5],
      row.names = NULL
    )
  })
  bind_rows(rows)
}

## ---- 3. Compute asymptotic diversity for both communities -----------------

cat("\nEstimating asymptotic diversity for", nrow(insect), "insect samples...\n")
set.seed(1)
insect_asymp <- compute_asymptotic_metrics(insect) %>% left_join(meta, by = "sample_id")
cat("Estimating asymptotic diversity for", nrow(fungal), "fungal samples (this takes a few minutes)...\n")
# Re-seed before the fungal pass so its ChaoShannon/ChaoSimpson bootstrap SEs
# (B=200) are reproducible independent of how many RNG draws the insect pass
# above consumed -- otherwise a change to the insect table (e.g. the
# 2026-09-09 singleton-filter fix) silently shifts every fungal SE/CI even
# though the fungal point estimates are untouched.
set.seed(1)
fungal_asymp <- compute_asymptotic_metrics(fungal) %>% left_join(meta, by = "sample_id")

metrics <- c("richness", "shannon", "simpson")
metric_labels <- c(richness = "Species richness", shannon = "Shannon diversity (Hill q=1)",
                    simpson = "Simpson diversity (Hill q=2)")

cat("\nMedian sample coverage (fraction of true community captured): insect =",
    round(median(insect_asymp$coverage), 3), ", fungal =", round(median(fungal_asymp$coverage), 3), "\n")

write.csv(insect_asymp, file.path(out_data_dir, "insect_asymptotic_diversity.csv"), row.names = FALSE)
write.csv(fungal_asymp, file.path(out_data_dir, "fungal_asymptotic_diversity.csv"), row.names = FALSE)

## ---- 4. Cross-community comparison: does insect diversity track fungal diversity? ----
# Spearman correlation per metric, both on OBSERVED and ASYMPTOTIC values, so
# the effect of extrapolation on the conclusion itself is visible directly --
# same samples used for both value types (fungal already restricted to the
# matched set in section 1, no rarefaction-survivor dropout as in the
# rarefy-and-average script).

combined <- insect_asymp %>%
  select(sample_id, site, lure, date, trap_id,
         insect_richness_observed = richness_observed, insect_richness_asymptotic = richness_asymptotic,
         insect_shannon_observed = shannon_observed, insect_shannon_asymptotic = shannon_asymptotic,
         insect_simpson_observed = simpson_observed, insect_simpson_asymptotic = simpson_asymptotic) %>%
  inner_join(
    fungal_asymp %>% select(sample_id,
                             fungal_richness_observed = richness_observed, fungal_richness_asymptotic = richness_asymptotic,
                             fungal_shannon_observed = shannon_observed, fungal_shannon_asymptotic = shannon_asymptotic,
                             fungal_simpson_observed = simpson_observed, fungal_simpson_asymptotic = simpson_asymptotic),
    by = "sample_id"
  )
cat("\n", nrow(combined), "samples with both insect and fungal diversity for cross-community comparison.\n")

value_types <- c("observed", "asymptotic")
cross_corr <- bind_rows(lapply(value_types, function(vt) {
  bind_rows(lapply(metrics, function(m) {
    x <- combined[[paste0("insect_", m, "_", vt)]]
    y <- combined[[paste0("fungal_", m, "_", vt)]]
    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    data.frame(value_type = vt, metric = m, rho = unname(ct$estimate), S_stat = unname(ct$statistic),
               p_value = ct$p.value, n = nrow(combined))
  }))
}))
cat("\n--- Insect vs. fungal diversity, per-sample Spearman correlation (observed and asymptotic) ---\n")
print(cross_corr, row.names = FALSE)
write.csv(cross_corr, file.path(out_data_dir, "asymptotic_diversity_cross_community_correlation.csv"), row.names = FALSE)

## ---- 4b. Does lure modify the insect~fungal relationship? (asymptotic values) ----
# Same test as insect_fungal_alpha_diversity.full_insect_table.r's section
# 5b (fungal ~ site+date+insect*lure, permutation F-test on the interaction
# term, lure permuted among traps within site), applied to the asymptotic
# estimates.

n_perm_interaction <- 999

test_interaction_lure <- function(fungal_y, insect_x, dat_base, trap_lure_map, n_perm = 999) {
  dat <- dat_base
  dat$y <- fungal_y
  dat$x <- insect_x
  fit_F <- function(d) {
    anova(lm(y ~ site + date + x + lure, data = d), lm(y ~ site + date + x * lure, data = d))[2, "F"]
  }
  obs_F <- fit_F(dat)
  perm_F <- replicate(n_perm, {
    perm_map <- trap_lure_map %>% group_by(site) %>% mutate(lure = sample(lure)) %>% ungroup()
    lure_lookup <- setNames(as.character(perm_map$lure), perm_map$trap_id)
    d2 <- dat; d2$lure <- lure_lookup[d2$trap_id]; fit_F(d2)
  })
  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)
  c(F = unname(obs_F), p_perm = p_perm)
}

trap_lure_map_combined <- combined %>% distinct(trap_id, site, lure)
cat("\nTesting insect_x:lure interaction on fungal asymptotic diversity (", n_perm_interaction, "permutations each)...\n")
interaction_tests <- bind_rows(lapply(metrics, function(m) {
  r <- test_interaction_lure(combined[[paste0("fungal_", m, "_asymptotic")]],
                              combined[[paste0("insect_", m, "_asymptotic")]],
                              combined, trap_lure_map_combined, n_perm_interaction)
  data.frame(metric = m, F = r["F"], p_perm = r["p_perm"], n = nrow(combined))
}))
rownames(interaction_tests) <- NULL
cat("\n--- Does lure modify the insect~fungal asymptotic-diversity relationship? ---\n")
print(interaction_tests, row.names = FALSE)
write.csv(interaction_tests, file.path(out_data_dir, "asymptotic_diversity_insect_fungal_interaction_by_lure.csv"),
          row.names = FALSE)

## ---- 5. Within-community site/lure/date effects on each metric (asymptotic values) ----
# Same permutation logic and linear+quadratic+cubic date test as
# insect_fungal_alpha_diversity.full_insect_table.r's section 6, applied to
# each community's own asymptotic diversity.

n_perm <- 999

test_term_date <- function(y, dat_base, n_perm = 999, block_var = "trap_id") {
  dat <- dat_base; dat$y <- y
  fit_model <- function(d) {
    d$date_c <- as.numeric(d$date) - mean(as.numeric(d$date))
    lm(y ~ site + lure + date_c + I(date_c^2) + I(date_c^3), data = d)
  }
  fit_stats <- function(d) {
    cs <- coef(summary(fit_model(d)))
    c(t_linear = cs["date_c", "t value"], t_quad = cs["I(date_c^2)", "t value"],
      t_cubic = cs["I(date_c^3)", "t value"])
  }
  obs <- fit_stats(dat)
  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)
  perm_stats <- t(apply(perm_ids, 1, function(idx) { d2 <- dat; d2$date <- dat$date[idx]; fit_stats(d2) }))
  p_linear <- (sum(abs(perm_stats[, "t_linear"]) >= abs(obs["t_linear"])) + 1) / (n_perm + 1)
  p_quad   <- (sum(abs(perm_stats[, "t_quad"]) >= abs(obs["t_quad"])) + 1) / (n_perm + 1)
  p_cubic  <- (sum(abs(perm_stats[, "t_cubic"]) >= abs(obs["t_cubic"])) + 1) / (n_perm + 1)

  full_aov <- anova(fit_model(dat))
  r2_of <- function(term) full_aov[term, "Sum Sq"] / (full_aov[term, "Sum Sq"] + full_aov["Residuals", "Sum Sq"])

  c(t_linear = unname(obs["t_linear"]), r2_linear = unname(r2_of("date_c")), p_linear = p_linear,
    t_quad = unname(obs["t_quad"]), r2_quad = unname(r2_of("I(date_c^2)")), p_quad = p_quad,
    t_cubic = unname(obs["t_cubic"]), r2_cubic = unname(r2_of("I(date_c^3)")), p_cubic = p_cubic)
}

test_term_lure <- function(y, dat_base, trap_lure_map, n_perm = 999) {
  dat <- dat_base; dat$y <- y
  fit_F <- function(d) anova(lm(y ~ site + date, data = d), lm(y ~ site + lure + date, data = d))[2, "F"]
  obs_F <- fit_F(dat)
  perm_F <- replicate(n_perm, {
    perm_map <- trap_lure_map %>% group_by(site) %>% mutate(lure = sample(lure)) %>% ungroup()
    lure_lookup <- setNames(as.character(perm_map$lure), perm_map$trap_id)
    d2 <- dat; d2$lure <- lure_lookup[d2$trap_id]; fit_F(d2)
  })
  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)
  full_aov <- anova(lm(y ~ site + lure + date, data = dat))
  r2 <- full_aov["lure", "Sum Sq"] / (full_aov["lure", "Sum Sq"] + full_aov["Residuals", "Sum Sq"])
  c(stat = unname(obs_F), r2 = unname(r2), p_perm = p_perm)
}

test_term_site <- function(y, dat_base, trap_site_map, n_perm = 999) {
  dat <- dat_base; dat$y <- y
  fit_F <- function(d) anova(lm(y ~ lure + date, data = d), lm(y ~ site + lure + date, data = d))[2, "F"]
  obs_F <- fit_F(dat)
  perm_F <- replicate(n_perm, {
    perm_map <- trap_site_map %>% mutate(site = sample(site))
    site_lookup <- setNames(as.character(perm_map$site), perm_map$trap_id)
    d2 <- dat; d2$site <- site_lookup[d2$trap_id]; fit_F(d2)
  })
  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)
  full_aov <- anova(lm(y ~ site + lure + date, data = dat))
  r2 <- full_aov["site", "Sum Sq"] / (full_aov["site", "Sum Sq"] + full_aov["Residuals", "Sum Sq"])
  c(stat = unname(obs_F), r2 = unname(r2), p_perm = p_perm)
}

communities <- list(insect = insect_asymp, fungal = fungal_asymp)

cat("\nTesting site/lure/date effects on each asymptotic metric (", n_perm, "permutations each)...\n")
covariate_tests <- bind_rows(lapply(names(communities), function(comm_name) {
  dat_base <- communities[[comm_name]]
  trap_lure_map <- dat_base %>% distinct(trap_id, site, lure)
  trap_site_map <- dat_base %>% distinct(trap_id, site)
  bind_rows(lapply(metrics, function(m) {
    y <- dat_base[[paste0(m, "_asymptotic")]]
    r_date <- test_term_date(y, dat_base, n_perm)
    r_lure <- test_term_lure(y, dat_base, trap_lure_map, n_perm)
    r_site <- test_term_site(y, dat_base, trap_site_map, n_perm)
    bind_rows(
      data.frame(community = comm_name, metric = m, term = "date_linear", stat_type = "t",
                 stat = r_date["t_linear"], r2 = r_date["r2_linear"], p_perm = r_date["p_linear"]),
      data.frame(community = comm_name, metric = m, term = "date_quadratic", stat_type = "t",
                 stat = r_date["t_quad"], r2 = r_date["r2_quad"], p_perm = r_date["p_quad"]),
      data.frame(community = comm_name, metric = m, term = "date_cubic", stat_type = "t",
                 stat = r_date["t_cubic"], r2 = r_date["r2_cubic"], p_perm = r_date["p_cubic"]),
      data.frame(community = comm_name, metric = m, term = "lure", stat_type = "F",
                 stat = r_lure["stat"], r2 = r_lure["r2"], p_perm = r_lure["p_perm"]),
      data.frame(community = comm_name, metric = m, term = "site", stat_type = "F",
                 stat = r_site["stat"], r2 = r_site["r2"], p_perm = r_site["p_perm"])
    )
  }))
}))
rownames(covariate_tests) <- NULL
cat("\n--- Site/lure/date effects on asymptotic diversity (each community's own table) ---\n")
print(covariate_tests, row.names = FALSE)
write.csv(covariate_tests, file.path(out_data_dir, "asymptotic_diversity_covariate_tests.csv"), row.names = FALSE)

date_shape <- covariate_tests %>%
  filter(term %in% c("date_linear", "date_quadratic", "date_cubic")) %>%
  select(community, metric, term, stat, p_perm) %>%
  pivot_wider(names_from = term, values_from = c(stat, p_perm)) %>%
  mutate(shape = case_when(
    p_perm_date_cubic < 0.10 ~ "complex (significant cubic term)",
    p_perm_date_quadratic < 0.10 & stat_date_quadratic < 0 ~ "hump (peaks mid-season)",
    p_perm_date_quadratic < 0.10 & stat_date_quadratic > 0 ~ "dip (troughs mid-season)",
    p_perm_date_linear < 0.10 & stat_date_linear > 0 ~ "linear increase (late-season)",
    p_perm_date_linear < 0.10 & stat_date_linear < 0 ~ "linear decrease (early-season)",
    TRUE ~ "no significant date pattern"
  ))
cat("\n--- Date shape classification, asymptotic diversity (linear/quadratic/cubic, cubic takes priority) ---\n")
print(date_shape, row.names = FALSE)
write.csv(date_shape, file.path(out_data_dir, "asymptotic_diversity_date_shape.csv"), row.names = FALSE)

## ---- 6. Figures --------------------------------------------------------------

alpha_theme <- theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))

# 6a. Observed vs. asymptotic (extrapolated) value by community -- the direct
# visualization of "how much more is out there" per metric.
long_asymp <- bind_rows(
  insect_asymp %>% mutate(community = "Insect"),
  fungal_asymp %>% mutate(community = "Fungal")
) %>%
  select(community, sample_id, all_of(paste0(rep(metrics, each = 2), c("_observed", "_asymptotic")))) %>%
  pivot_longer(-c(community, sample_id), names_to = c("metric", "value_type"),
               names_pattern = "(.*)_(observed|asymptotic)", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         value_type = factor(value_type, levels = c("observed", "asymptotic"),
                              labels = c("Observed", "Asymptotic (Chao-extrapolated)")))

dist_plot <- ggplot(long_asymp, aes(x = community, y = value, fill = value_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.1, shape = 21) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c("Observed" = "grey70", "Asymptotic (Chao-extrapolated)" = "#2166ac"), name = NULL) +
  labs(x = NULL, y = NULL, title = "Observed vs. asymptotic diversity by community (full insect table)") +
  alpha_theme
ggsave(file.path(out_fig_dir, "asymptotic_diversity_by_community.png"), dist_plot, width = 10, height = 4, dpi = 150)

# 6b. Paired observed -> asymptotic per sample (richness only, most direct read
# on "how many undetected taxa are estimated").
paired_df <- bind_rows(
  insect_asymp %>% mutate(community = "Insect"),
  fungal_asymp %>% mutate(community = "Fungal")
) %>% select(community, sample_id, coverage, richness_observed, richness_asymptotic) %>%
  pivot_longer(c(richness_observed, richness_asymptotic), names_to = "value_type", values_to = "richness") %>%
  mutate(value_type = factor(value_type, levels = c("richness_observed", "richness_asymptotic"),
                              labels = c("Observed", "Asymptotic")))
paired_plot <- ggplot(paired_df, aes(x = value_type, y = richness, group = sample_id, color = coverage)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 1.6) +
  facet_wrap(~community, scales = "free_y") +
  scale_color_viridis_c(name = "Sample\ncoverage") +
  labs(x = NULL, y = "Species richness",
       title = "Per-sample observed vs. asymptotic (Chao-extrapolated) richness") +
  alpha_theme
ggsave(file.path(out_fig_dir, "asymptotic_richness_observed_vs_extrapolated.png"), paired_plot,
       width = 8, height = 4.5, dpi = 150)

# 6c. Site/lure effects, asymptotic values
site_lure_df <- bind_rows(
  insect_asymp %>% mutate(community = "Insect"),
  fungal_asymp %>% mutate(community = "Fungal")
) %>%
  select(community, site, lure, all_of(paste0(metrics, "_asymptotic"))) %>%
  rename_with(~ sub("_asymptotic$", "", .x), all_of(paste0(metrics, "_asymptotic"))) %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         community = factor(community, levels = c("Insect", "Fungal")))

site_lure_plot <- ggplot(site_lure_df, aes(x = site, y = value, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.2) +
  facet_wrap(vars(metric, community), nrow = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2", name = "Lure") +
  labs(x = "Site", y = NULL, title = "Asymptotic diversity by site and lure (full insect table)") +
  alpha_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_fig_dir, "asymptotic_diversity_by_site_lure.png"), site_lure_plot, width = 9, height = 8, dpi = 150)

# 6d. Date effect, asymptotic values
date_df <- bind_rows(
  insect_asymp %>% mutate(community = "Insect"),
  fungal_asymp %>% mutate(community = "Fungal")
) %>%
  select(community, site, date, all_of(paste0(metrics, "_asymptotic"))) %>%
  rename_with(~ sub("_asymptotic$", "", .x), all_of(paste0(metrics, "_asymptotic"))) %>%
  pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         community = factor(community, levels = c("Insect", "Fungal")))

date_plot <- ggplot(date_df, aes(x = date, y = value, fill = site)) +
  geom_point(shape = 21, size = 2, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.5) +
  facet_wrap(vars(metric, community), nrow = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", name = "Site") +
  labs(x = "Collection date", y = NULL, title = "Asymptotic diversity vs. collection date (full insect table)") +
  alpha_theme
ggsave(file.path(out_fig_dir, "asymptotic_diversity_by_date.png"), date_plot, width = 9, height = 8, dpi = 150)

# 6e. Cross-community correlation scatter, asymptotic values, one panel per metric
make_corr_panel <- function(m) {
  df <- combined
  x_col <- paste0("insect_", m, "_asymptotic"); y_col <- paste0("fungal_", m, "_asymptotic")
  cc <- cross_corr %>% filter(metric == m, value_type == "asymptotic")
  label <- paste0("rho = ", round(cc$rho, 2), "\np = ", signif(cc$p_value, 2))
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = as.numeric(date), shape = site)) +
    geom_point(size = 2.6, color = "black", stroke = 0.3) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                          midpoint = median(as.numeric(df$date)), name = "Date") +
    scale_shape_manual(values = c(21, 22, 23, 24), name = "Site") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = label, size = 3.2,
             fill = alpha("white", 0.75)) +
    labs(x = paste0("Insect asymptotic ", metric_labels[m]), y = paste0("Fungal asymptotic ", metric_labels[m])) +
    alpha_theme
}
corr_panels <- lapply(metrics, make_corr_panel)
corr_grid <- gridExtra::grid.arrange(grobs = corr_panels, ncol = 3,
                                      top = "Insect (full table) vs. fungal asymptotic diversity (per-sample)")
ggsave(file.path(out_fig_dir, "asymptotic_diversity_cross_community_correlation.png"), corr_grid,
       width = 13, height = 4.5, dpi = 150)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
