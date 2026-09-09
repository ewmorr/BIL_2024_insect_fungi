##############################################################################
# Insect vs. Fungal Alpha Diversity -- Full Insect Table Variant
#
# Same analysis as insect_fungal_alpha_diversity.r (kept as-is for
# reference), except the insect community table is NOT restricted to
# Curculionidae + Latridiidae (the two families used throughout the rest of
# this project's insect/fungal comparisons). Here the insect table is every
# taxon in insect_species_tab.csv (98 families, 426 taxa total), all families.
#
# Insect-table filtering (changed 2026-09-09): RAW counts, empty-sample drop
# only -- NO global singleton filter (colSums > 1). This script previously
# applied that filter as noise-floor cleanup (the convention documented in
# project_organization.md's "When to prevalence-filter vs. not" for
# alpha-diversity work). It is dropped here for consistency with
# insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r, which MUST use
# raw counts (Chao-type estimators are built on the observed singleton count
# f1 -- see that script's header). With the observed script also on the raw
# table, "observed richness" now means the same thing in both, so the
# observed-vs-asymptotic comparison in that script is like-for-like. Effect
# is modest: observed richness rises a median of ~2 taxa/sample (Shannon and
# Simpson dominance barely move -- global singletons carry negligible p_i).
# See the 2026-09-09 entry in iterative_analysis_updates.md. (Lineage A's
# family-filtered insect_fungal_alpha_diversity.r is frozen as the backup and
# is NOT changed -- it keeps the singleton filter.)
#
# All prior comparisons in this project (CoCA, per-taxon date/lure tests,
# Procrustes) ask a beta-diversity question -- does community COMPOSITION
# covary. This script asks the alpha-diversity question instead: do the two
# communities' per-sample Shannon diversity, Simpson dominance, and richness
# (observed taxon count, no extrapolation) covary, and does each community's
# own diversity respond to site/lure/date the way its composition does?
#
# Community tables:
# - Insect: the full raw trap-catch table (all families, raw counts, only
#   all-zero sample rows dropped -- NO singleton filter, see the filtering
#   note above). No rarefaction -- trap catch counts aren't subject to the
#   same sequencing-depth artifact as amplicon reads.
# - Fungal: unchanged from insect_fungal_alpha_diversity.r -- the same
#   multiple-subsampling (rarefy to depth 5000 x 100 iterations) approach
#   used for the fungal NMDS/PERMANOVA (fungal_community_seasonality.r,
#   fig2_nmds_procrustes.R), averaging per-sample alpha diversity across
#   iterations instead of averaging Bray-Curtis distance matrices.
#
# Two comparisons, per user request:
# 1. Per-sample insect-vs-fungal correlation for each metric (Spearman) --
#    the alpha-diversity analog of the whole-community Procrustes check.
# 2. Site/lure/date effects on each metric within each community, using the
#    same permutation schemes already established in this project: date
#    permuted within trap (fungal_community_seasonality.r), lure permuted
#    among traps within site (insect_pcoa_lure_association.r), and site
#    permuted freely among all traps (same logic as the lure test, one level
#    up -- site, like lure, is a trap-level constant, so it is tested by
#    permuting the trap-level label rather than individual samples).
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(gridExtra)
source("library/library.R")

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_alpha_diversity_full_insect_table"
out_fig_dir <- "figures/compare_insects_fungi_alpha_diversity_full_insect_table"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data (full insect table, all families) ----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

cat("Full insect table:", nrow(sp_tab), "taxa across", length(unique(sp_tab$Family)), "families.\n")
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

# NO singleton filter (changed 2026-09-09, see header) -- raw counts, drop
# only any all-zero sample rows.
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

# Restrict to ASVs confirmed as Fungi (same filter as every other script here).
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
    ncol(fungal), "fungal ASVs.\n")

## ---- 2. Alpha-diversity helper: Shannon, Simpson dominance, richness -------
# Simpson DOMINANCE (D = sum(p_i^2), higher = more dominated by few taxa) is
# computed directly rather than via vegan's index="simpson" (which returns
# 1-D, a diversity rather than a dominance measure) so the sign matches what
# was requested. Richness is the observed taxon count (specnumber) -- no
# extrapolated/Chao-type estimator.

compute_alpha_metrics <- function(x) {
  x <- as.matrix(x)
  props <- sweep(x, 1, rowSums(x), "/")
  data.frame(
    sample_id = rownames(x),
    shannon = vegan::diversity(x, index = "shannon"),
    simpson_dominance = rowSums(props^2),
    richness = vegan::specnumber(x),
    row.names = NULL
  )
}

## ---- 3. Insect alpha diversity (full raw trap-catch table) -----------------

insect_alpha <- compute_alpha_metrics(insect) %>% left_join(meta, by = "sample_id")
cat("\nInsect alpha diversity computed for", nrow(insect_alpha), "samples.\n")

## ---- 4. Fungal alpha diversity (rarefy-and-average, same as the ordination) ----

rarefy_depth <- 5000
rarefy_iterations <- 100

cat("\nRarefying fungal ASV table (", ncol(fungal), "ASVs,", nrow(fungal),
    "samples) to depth", rarefy_depth, "x", rarefy_iterations, "iterations...\n")
fungal_rare_list <- multiple_subsamples(x = fungal, depth = rarefy_depth, iterations = rarefy_iterations)
n_dropped <- nrow(fungal) - nrow(fungal_rare_list[[1]])
if (n_dropped > 0) {
  dropped_ids <- setdiff(rownames(fungal), rownames(fungal_rare_list[[1]]))
  cat(n_dropped, "sample(s) dropped (total reads fell below depth", rarefy_depth, "):",
      paste(dropped_ids, collapse = ", "), "\n")
}

fungal_alpha_iterations <- bind_rows(lapply(fungal_rare_list, compute_alpha_metrics), .id = "iteration")

# Average each metric across the 100 rarefaction iterations per sample -- the
# alpha-diversity analog of avg_matrix_list() averaging distance matrices.
fungal_alpha <- fungal_alpha_iterations %>%
  group_by(sample_id) %>%
  summarize(shannon = mean(shannon), simpson_dominance = mean(simpson_dominance),
            richness = mean(richness), .groups = "drop") %>%
  left_join(meta, by = "sample_id")
cat("Fungal alpha diversity computed for", nrow(fungal_alpha), "samples",
    "(averaged across", rarefy_iterations, "rarefaction iterations).\n")

metrics <- c("shannon", "simpson_dominance", "richness")
metric_labels <- c(shannon = "Shannon diversity", simpson_dominance = "Simpson dominance",
                    richness = "Richness (observed taxa)")

write.csv(insect_alpha, file.path(out_data_dir, "insect_alpha_diversity.csv"), row.names = FALSE)
write.csv(fungal_alpha, file.path(out_data_dir, "fungal_alpha_diversity.csv"), row.names = FALSE)

## ---- 5. Cross-community comparison: does insect diversity track fungal diversity? ----
# Spearman correlation per metric, on the samples both communities share
# (fungal rarefaction survivors) -- the alpha-diversity analog of the
# whole-community Procrustes check (insect_fungal_procrustes.r), which also
# used a plain (unconstrained) permutation/correlation test rather than
# trap-blocking.

combined <- insect_alpha %>%
  select(sample_id, site, lure, date, trap_id, insect_shannon = shannon,
         insect_simpson_dominance = simpson_dominance, insect_richness = richness) %>%
  inner_join(
    fungal_alpha %>% select(sample_id, fungal_shannon = shannon,
                             fungal_simpson_dominance = simpson_dominance, fungal_richness = richness),
    by = "sample_id"
  )
cat("\n", nrow(combined), "samples with both insect and fungal alpha diversity for cross-community comparison.\n")

cross_corr <- bind_rows(lapply(metrics, function(m) {
  ct <- suppressWarnings(cor.test(combined[[paste0("insect_", m)]], combined[[paste0("fungal_", m)]],
                                   method = "spearman"))
  data.frame(metric = m, rho = unname(ct$estimate), S_stat = unname(ct$statistic),
             p_value = ct$p.value, n = nrow(combined))
}))
cat("\n--- Insect vs. fungal alpha diversity, per-sample Spearman correlation ---\n")
print(cross_corr, row.names = FALSE)
write.csv(cross_corr, file.path(out_data_dir, "alpha_diversity_cross_community_correlation.csv"), row.names = FALSE)

## ---- 5b. Does lure modify the insect~fungal relationship? ------------------
# The pooled Spearman correlation above (and fig10_insect_fungal_alpha_
# diversity_correlation.R, which visualizes it split by lure with a purely
# descriptive per-lure Spearman rho/p in each panel) doesn't test whether
# lure actually changes that relationship -- splitting the data 3 ways for a
# plot isn't the same as testing an interaction. This tests it directly:
# fungal_y ~ site + date + insect_x * lure, i.e. does the insect_x:lure
# interaction term explain more variance than a model with only the insect_x
# and lure main effects (plus site/date controls, matching every other test
# in this script)? Tested via the same trap-level permutation this project
# uses for lure everywhere else (test_term_lure below, insect_pcoa_lure_
# association.r): lure permuted among the traps WITHIN each site, since lure
# is constant within trap (one trap per site x lure combination).
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
cat("\nTesting insect_x:lure interaction on fungal diversity (", n_perm_interaction, "permutations each)...\n")
interaction_tests <- bind_rows(lapply(metrics, function(m) {
  r <- test_interaction_lure(combined[[paste0("fungal_", m)]], combined[[paste0("insect_", m)]],
                              combined, trap_lure_map_combined, n_perm_interaction)
  data.frame(metric = m, F = r["F"], p_perm = r["p_perm"], n = nrow(combined))
}))
rownames(interaction_tests) <- NULL
cat("\n--- Does lure modify the insect~fungal relationship? (fungal ~ site+date+insect*lure) ---\n")
print(interaction_tests, row.names = FALSE)
write.csv(interaction_tests, file.path(out_data_dir, "alpha_diversity_insect_fungal_interaction_by_lure.csv"),
          row.names = FALSE)

## ---- 6. Within-community site/lure/date effects on each metric -------------
# Same permutation logic used throughout this project for testing trap-level
# fixed effects, applied here to a univariate diversity metric via lm()
# instead of a per-taxon abundance or a distance matrix:
#   - date: varies within trap (repeated sampling over the season) -- permute
#     date within trap, holding site/lure fixed (fungal_community_seasonality.r).
#   - lure: constant within trap, one trap per site x lure combination --
#     permute lure among the 3 traps WITHIN each site (insect_pcoa_lure_
#     association.r); this is the permutation that matches how lure actually
#     varies in the design.
#   - site: also constant within trap, but with no grouping level above it --
#     tested the same way as lure, one level up: permute the site label
#     freely among ALL traps.
# As with the community-level PERMANOVA tables elsewhere in this project,
# these per-term p-values are reported without multiple-testing correction
# (they're a small, fixed set of a priori terms, not a taxon screen).
#
# Date is tested with linear + quadratic + cubic terms in ONE model (date
# centered to avoid a huge-magnitude collinear polynomial term), extending
# this project's established linear+quadratic date pattern (fungal_
# community_seasonality.r, insect_fungal_coca_analysis.general_workflow.r --
# see project_organization.md's "Linear+quadratic date test pattern") one
# order further to cubic. Added 2026-09-01 specifically for this alpha-
# diversity check: fig9_alpha_diversity_by_lure.R's loess trend lines showed
# several panels (e.g. insect Shannon diversity under Ethanol, several
# fungal Simpson-dominance panels) with more structure than a single hump/
# dip can capture, so a cubic term is tested here as a targeted follow-up --
# this is NOT (yet) a project-wide convention the way linear+quadratic is.

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

communities <- list(insect = insect_alpha, fungal = fungal_alpha)

cat("\nTesting site/lure/date effects on each metric (", n_perm, "permutations each)...\n")
covariate_tests <- bind_rows(lapply(names(communities), function(comm_name) {
  dat_base <- communities[[comm_name]]
  trap_lure_map <- dat_base %>% distinct(trap_id, site, lure)
  trap_site_map <- dat_base %>% distinct(trap_id, site)
  bind_rows(lapply(metrics, function(m) {
    y <- dat_base[[m]]
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
cat("\n--- Site/lure/date effects on alpha diversity (each community's own table) ---\n")
print(covariate_tests, row.names = FALSE)
write.csv(covariate_tests, file.path(out_data_dir, "alpha_diversity_covariate_tests.csv"), row.names = FALSE)

# Shape classification per community x metric, same cascading logic as
# fungal_community_seasonality.r's quadratic-takes-priority-over-linear rule
# (see project_organization.md), extended one level further: a significant
# cubic term takes priority over quadratic, which takes priority over linear.
# p_perm (uncorrected) is used directly, same as the rest of this covariate
# table -- these are 2 communities x 3 metrics = 6 date tests total, the
# same kind of small fixed a priori set as the site/lure/date PERMANOVA
# terms elsewhere in this project, not a taxon screen needing BH-FDR.
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
cat("\n--- Date shape classification (linear/quadratic/cubic, cubic takes priority) ---\n")
print(date_shape, row.names = FALSE)
write.csv(date_shape, file.path(out_data_dir, "alpha_diversity_date_shape.csv"), row.names = FALSE)

## ---- 7. Figures --------------------------------------------------------------

alpha_theme <- theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))

# 7a. Distribution of each metric by community (descriptive)
long_alpha <- bind_rows(
  insect_alpha %>% mutate(community = "Insect") %>% select(community, sample_id, all_of(metrics)),
  fungal_alpha %>% mutate(community = "Fungal") %>% select(community, sample_id, all_of(metrics))
) %>% pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels))

dist_plot <- ggplot(long_alpha, aes(x = community, y = value, fill = community)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.6, shape = 21, color = "black") +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(Insect = "#2166ac", Fungal = "#b2182b"), guide = "none") +
  labs(x = NULL, y = NULL, title = "Alpha diversity by community (full insect table)") +
  alpha_theme
ggsave(file.path(out_fig_dir, "alpha_diversity_by_community.png"), dist_plot, width = 9, height = 4, dpi = 150)

# 7b. Site/lure effects: boxplot by site, colored by lure, faceted by community x metric
site_lure_df <- bind_rows(
  insect_alpha %>% mutate(community = "Insect"),
  fungal_alpha %>% mutate(community = "Fungal")
) %>% pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         community = factor(community, levels = c("Insect", "Fungal")))

# facet_wrap (not facet_grid) so each metric x community panel gets its own
# y-scale -- insect and fungal richness are on very different scales, and
# facet_grid(scales="free_y") only frees scales across rows, not across
# columns within a row.
site_lure_plot <- ggplot(site_lure_df, aes(x = site, y = value, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.2) +
  facet_wrap(vars(metric, community), nrow = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2", name = "Lure") +
  labs(x = "Site", y = NULL, title = "Alpha diversity by site and lure (full insect table)") +
  alpha_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_fig_dir, "alpha_diversity_by_site_lure.png"), site_lure_plot, width = 9, height = 8, dpi = 150)

# 7c. Date effect: scatter of value vs. date, colored by site, faceted by community x metric
date_plot <- ggplot(site_lure_df, aes(x = date, y = value, fill = site)) +
  geom_point(shape = 21, size = 2, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.5) +
  facet_wrap(vars(metric, community), nrow = 3, scales = "free_y") +
  scale_fill_brewer(palette = "Set2", name = "Site") +
  labs(x = "Collection date", y = NULL, title = "Alpha diversity vs. collection date (full insect table)") +
  alpha_theme
ggsave(file.path(out_fig_dir, "alpha_diversity_by_date.png"), date_plot, width = 9, height = 8, dpi = 150)

# 7d. Cross-community correlation scatter, one panel per metric
make_corr_panel <- function(m) {
  df <- combined
  x_col <- paste0("insect_", m); y_col <- paste0("fungal_", m)
  cc <- cross_corr %>% filter(metric == m)
  label <- paste0("rho = ", round(cc$rho, 2), "\np = ", signif(cc$p_value, 2))
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = as.numeric(date), shape = site)) +
    geom_point(size = 2.6, color = "black", stroke = 0.3) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                          midpoint = median(as.numeric(df$date)), name = "Date") +
    scale_shape_manual(values = c(21, 22, 23, 24), name = "Site") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = label, size = 3.2,
             fill = alpha("white", 0.75)) +
    labs(x = paste0("Insect ", metric_labels[m]), y = paste0("Fungal ", metric_labels[m])) +
    alpha_theme
}
corr_panels <- lapply(metrics, make_corr_panel)
corr_grid <- gridExtra::grid.arrange(grobs = corr_panels, ncol = 3,
                                      top = "Insect (full table) vs. fungal alpha diversity (per-sample)")
ggsave(file.path(out_fig_dir, "alpha_diversity_cross_community_correlation.png"), corr_grid,
       width = 13, height = 4.5, dpi = 150)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
