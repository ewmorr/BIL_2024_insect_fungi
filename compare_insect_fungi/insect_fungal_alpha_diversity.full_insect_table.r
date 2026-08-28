##############################################################################
# Insect vs. Fungal Alpha Diversity -- Full Insect Table Variant
#
# Same analysis as insect_fungal_alpha_diversity.r (kept as-is for
# reference), except the insect community table is NOT restricted to
# Curculionidae + Latridiidae (the two families used throughout the rest of
# this project's insect/fungal comparisons). Here the insect table is every
# taxon in insect_species_tab.csv (79 families, 426 taxa total) -- the same
# singleton-taxon and empty-sample filters are applied, just without the
# family restriction.
#
# All prior comparisons in this project (CoCA, per-taxon date/lure tests,
# Procrustes) ask a beta-diversity question -- does community COMPOSITION
# covary. This script asks the alpha-diversity question instead: do the two
# communities' per-sample Shannon diversity, Simpson dominance, and richness
# (observed taxon count, no extrapolation) covary, and does each community's
# own diversity respond to site/lure/date the way its composition does?
#
# Community tables:
# - Insect: the full raw trap-catch table (all families, singleton taxa and
#   resulting empty samples dropped -- same filters insect_ords.R/fig2 apply
#   to the Curculionidae+Latridiidae subset, just not restricted to those
#   two families first). No rarefaction -- trap catch counts aren't subject
#   to the same sequencing-depth artifact as amplicon reads.
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

# drop singleton taxa (present as a single individual total) and any
# resulting all-zero sample rows -- same filters applied to the
# Curculionidae+Latridiidae subset elsewhere in this project.
sp_tab.t[, colSums(sp_tab.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]
cat("After dropping singleton taxa and empty samples:", nrow(insect_full), "samples,",
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

n_perm <- 999

test_term_date <- function(y, dat_base, n_perm = 999, block_var = "trap_id") {
  dat <- dat_base; dat$y <- y
  fit_stat <- function(d) coef(summary(lm(y ~ site + lure + date, data = d)))["date", "t value"]
  obs_t <- fit_stat(dat)
  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)
  perm_t <- apply(perm_ids, 1, function(idx) { d2 <- dat; d2$date <- dat$date[idx]; fit_stat(d2) })
  p_perm <- (sum(abs(perm_t) >= abs(obs_t)) + 1) / (n_perm + 1)
  full_aov <- anova(lm(y ~ site + lure + date, data = dat))
  r2 <- full_aov["date", "Sum Sq"] / (full_aov["date", "Sum Sq"] + full_aov["Residuals", "Sum Sq"])
  c(stat = unname(obs_t), r2 = unname(r2), p_perm = p_perm)
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
      data.frame(community = comm_name, metric = m, term = "date", stat_type = "t",
                 stat = r_date["stat"], r2 = r_date["r2"], p_perm = r_date["p_perm"]),
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
