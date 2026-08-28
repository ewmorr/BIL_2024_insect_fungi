##############################################################################
# Insect Community PCoA1-3 -- Association with Date, Site, and Lure
# (Prevalence-Filtered Insect Table Variant)
#
# insect_fungal_coca_analysis.prevalence_filtered_insect.r reran the CoCA
# workflow on the >=5-sample individual-taxon prevalence-filtered insect
# table (153 taxa, all families -- insect_exploratory/
# insect_ords.prevalence_filter.R) instead of the Curculionidae+Latridiidae
# restriction, and found a new insect_PCoA2 <-> fungal association (6/200
# candidates significant without date, 7/200 with date, same taxa/sign
# either way) that was completely null (0/200) under the family-filtered
# table throughout every earlier check in this project. This script asks two
# follow-up questions:
#
# Part 1: what IS insect_PCoA2 under the broader table -- does it track
#   date, site, or lure? Mirrors insect_pcoa_date_lure_association.r exactly
#   (same two permutation schemes: date permuted within trap, lure permuted
#   among traps within site; site gets a parametric-only test), just rebuilt
#   on the prevalence-filtered insect table.
# Part 2: which fungal taxa are driving the insect_PCoA2 CoCA hit, and what
#   are they taxonomically? Pulls the significant-hit taxon lists straight
#   from insect_fungal_coca_analysis.prevalence_filtered_insect.r's saved
#   output and joins to ASVs_taxonomy.tsv.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("vegan", "dplyr", "tidyr", "tibble", "ggplot2")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"
out_fig_dir <- "figures/compare_insects_fungi_top3axes_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data (prevalence-filtered insect table) ----
# Fungal data is loaded only to restrict to the same sample set used in the
# prevalence-filtered CoCA script -- fungal read counts themselves are never
# used below.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]  # one all-zero row with no ID
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]
cat(sum(keep_taxa), "of", ncol(sp_tab.t), "insect taxa retained at >=5-sample prevalence (all families).\n")

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full), "insect samples have matching fungal ASV data.\n")
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (>=5-sample prevalence filter, all families).\n")

# Design checks the two permutation schemes below depend on.
trap_lure_map <- meta %>% distinct(trap_id, site, lure)
stopifnot(all(table(trap_lure_map$site, trap_lure_map$lure) == 1))
cat(nrow(trap_lure_map), "traps, one per site x lure combination.\n")

dates_per_trap <- meta %>% distinct(trap_id, date) %>% count(trap_id) %>% pull(n)
cat("Dates per trap: min", min(dates_per_trap), "max", max(dates_per_trap),
    "-- within-trap date permutation assumes every trap was visited on (close to) the same date set.\n")

## ---- 2. Insect PCoA1-3 (same construction as the CoCA script) -------------
# Hellinger transform + Euclidean distance = Hellinger distance -- see
# insect_fungal_coca_analysis.general_workflow.r for the reasoning.

insect_hel <- decostand(insect, method = "hellinger")
insect_dist <- vegdist(insect_hel, method = "euclidean")
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)

pos_eig <- insect_pcoa$eig[insect_pcoa$eig > 0]
var_explained <- pos_eig / sum(pos_eig)
cat("Variance explained, PCoA1-3:", paste(round(100 * head(var_explained, 3), 1), collapse = "%, "), "%\n")

insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(insect)))

## ---- 3. Per-axis date and lure association tests --------------------------
# Identical logic/formulas to insect_pcoa_date_lure_association.r -- see that
# script's header for the full rationale behind each permutation scheme and
# the partial/semi-partial R2 definitions.

n_perm <- 999

term_r2s <- function(dat, term) {
  rhs <- c(setdiff(c("site", "lure", "date"), term), term)
  full_aov <- anova(lm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = dat))
  ss_term <- full_aov[term, "Sum Sq"]
  ss_resid <- full_aov["Residuals", "Sum Sq"]
  ss_total <- sum((dat$y - mean(dat$y))^2)
  c(r2_partial = ss_term / (ss_term + ss_resid), r2_semipartial = ss_term / ss_total)
}

test_axis_date <- function(axis_vals, dat_base, n_perm) {
  dat <- dat_base
  dat$y <- axis_vals

  fit_F <- function(d) {
    m0 <- lm(y ~ site + lure, data = d)
    m1 <- lm(y ~ site + lure + date, data = d)
    anova(m0, m1)[2, "F"]
  }

  obs_F <- fit_F(dat)

  perm_F <- replicate(n_perm, {
    d2 <- dat %>% group_by(trap_id) %>% mutate(date = sample(date)) %>% ungroup()
    fit_F(d2)
  })

  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)

  c(F_stat = obs_F, term_r2s(dat, "date"),
    r_pearson = cor(dat$y, as.numeric(dat$date)), p_perm = p_perm)
}

test_axis_lure <- function(axis_vals, dat_base, trap_lure_map, n_perm) {
  dat <- dat_base
  dat$y <- axis_vals

  fit_F <- function(d) {
    m0 <- lm(y ~ site + date, data = d)
    m1 <- lm(y ~ site + lure + date, data = d)
    anova(m0, m1)[2, "F"]
  }

  obs_F <- fit_F(dat)

  perm_F <- replicate(n_perm, {
    perm_map <- trap_lure_map %>%
      group_by(site) %>%
      mutate(lure = sample(lure)) %>%
      ungroup()
    lure_lookup <- setNames(as.character(perm_map$lure), perm_map$trap_id)
    d2 <- dat
    d2$lure <- lure_lookup[d2$trap_id]
    fit_F(d2)
  })

  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)

  c(F_stat = obs_F, term_r2s(dat, "lure"), p_perm = p_perm)
}

test_axis_site <- function(axis_vals, dat_base) {
  dat <- dat_base
  dat$y <- axis_vals
  m0 <- lm(y ~ lure + date, data = dat)
  m1 <- lm(y ~ site + lure + date, data = dat)
  a <- anova(m0, m1)
  c(F_stat = a[2, "F"], term_r2s(dat, "site"), p_param = a[2, "Pr(>F)"])
}

axis_names <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")

cat("\nTesting date and lure association for", paste(axis_names, collapse = ", "),
    "(", n_perm, "permutations each)...\n")

results <- bind_rows(lapply(axis_names, function(ax) {
  axis_num <- as.integer(sub("insect_PCoA", "", ax))
  y <- dat_base[[ax]]

  date_res <- test_axis_date(y, dat_base, n_perm)
  lure_res <- test_axis_lure(y, dat_base, trap_lure_map, n_perm)
  site_res <- test_axis_site(y, dat_base)

  bind_rows(
    data.frame(axis = ax, pct_var_explained = 100 * var_explained[axis_num], term = "date",
               F_stat = date_res["F_stat"], r2_partial = date_res["r2_partial"],
               r2_semipartial = date_res["r2_semipartial"],
               r_pearson = date_res["r_pearson"], p_value = date_res["p_perm"], p_type = "permutation"),
    data.frame(axis = ax, pct_var_explained = 100 * var_explained[axis_num], term = "lure",
               F_stat = lure_res["F_stat"], r2_partial = lure_res["r2_partial"],
               r2_semipartial = lure_res["r2_semipartial"],
               r_pearson = NA, p_value = lure_res["p_perm"], p_type = "permutation"),
    data.frame(axis = ax, pct_var_explained = 100 * var_explained[axis_num], term = "site",
               F_stat = site_res["F_stat"], r2_partial = site_res["r2_partial"],
               r2_semipartial = site_res["r2_semipartial"],
               r_pearson = NA, p_value = site_res["p_param"], p_type = "parametric")
  )
}))
rownames(results) <- NULL

results <- results %>%
  group_by(term) %>%
  mutate(q_value = if (p_type[1] == "permutation") p.adjust(p_value, method = "BH") else NA_real_) %>%
  ungroup() %>%
  arrange(axis, term)

cat("\n--- insect_PCoA1-3 (prevalence-filtered table) vs. date / site / lure ---\n")
print(as.data.frame(results), row.names = FALSE)
write.csv(results, file.path(out_data_dir, "insect_pcoa_axes.date_site_lure_association.csv"), row.names = FALSE)

## ---- 4. Diagnostic plots ---------------------------------------------------

plot_dat <- dat_base %>% select(all_of(axis_names), site, lure, date) %>%
  pivot_longer(all_of(axis_names), names_to = "axis", values_to = "value")

date_plot <- ggplot(plot_dat, aes(x = date, y = value, color = site)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  facet_wrap(~axis, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Collection date", y = "Axis value",
       title = "insect_PCoA1-3 vs. collection date (prevalence-filtered insect table)") +
  theme_bw()
ggsave(file.path(out_fig_dir, "insect_pcoa_top3axes.by_date.png"), date_plot, width = 10, height = 4, dpi = 150)

lure_plot <- ggplot(plot_dat, aes(x = lure, y = value, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, aes(shape = site), size = 2) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(~axis, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Lure", y = "Axis value", title = "insect_PCoA1-3 by lure (prevalence-filtered insect table)") +
  theme_bw()
ggsave(file.path(out_fig_dir, "insect_pcoa_top3axes.by_lure.png"), lure_plot, width = 10, height = 4, dpi = 150)

## ---- 5. Taxonomy of the fungal taxa driving the insect_PCoA2 CoCA hit -----
# Pulls the significant-hit taxon lists straight from
# insect_fungal_coca_analysis.prevalence_filtered_insect.r's saved output
# (same candidate set, same per-taxon permutation test) and joins to
# ASVs_taxonomy.tsv, same taxonomy-prefix-stripping convention as
# presentation_items/fig4_association_taxonomic_breakdown.R.

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE) %>%
  rownames_to_column("taxon") %>%
  mutate(across(Kingdom:Species, ~ sub("^[a-z]__", "", .)))

pcoa2_only <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA2_only.csv"))
pcoa2_date <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA2_plus_date.csv"))

annotate_hits <- function(res, label) {
  res %>%
    filter(q_value < 0.10) %>%
    left_join(fungal_taxonomy, by = "taxon") %>%
    transmute(taxon, t_stat, q_value, coca_strength, Phylum, Class, Order, Family, Genus, Species,
              model = label) %>%
    arrange(q_value)
}

pcoa2_hits_no_date <- annotate_hits(pcoa2_only, "PCoA2_only (no date)")
pcoa2_hits_with_date <- annotate_hits(pcoa2_date, "PCoA2_plus_date")

cat("\n--- Taxonomy of fungal taxa significantly associated with insect_PCoA2 (no date), q<0.10 ---\n")
print(pcoa2_hits_no_date, row.names = FALSE)

cat("\n--- Taxonomy of fungal taxa significantly associated with insect_PCoA2 (date-adjusted), q<0.10 ---\n")
print(pcoa2_hits_with_date, row.names = FALSE)

shared_taxa <- intersect(pcoa2_hits_no_date$taxon, pcoa2_hits_with_date$taxon)
cat("\n", length(shared_taxa), "of", nrow(pcoa2_hits_no_date), "(no-date) /", nrow(pcoa2_hits_with_date),
    "(date-adjusted) hits are the same taxa in both models:", paste(shared_taxa, collapse = ", "), "\n")

pcoa2_taxonomy_summary <- bind_rows(pcoa2_hits_no_date, pcoa2_hits_with_date)
write.csv(pcoa2_taxonomy_summary, file.path(out_data_dir, "PCoA2_significant_taxa_taxonomy.csv"), row.names = FALSE)

cat("\nFamily tally across both PCoA2 hit lists (taxa counted once per model they appear in):\n")
print(pcoa2_taxonomy_summary %>% count(Family, sort = TRUE))

## ---- 6. Non-linear (quadratic) date association ---------------------------
# The insect_PCoA1-3 vs. date scatterplot (Part 4) shows insect_PCoA2 rising
# from early May to a peak around late May/early June, then falling back
# through July -- a hump shape, not a monotonic trend. The linear test in
# Part 3 (r=0.06, p=0.465) has essentially no power to detect that: a
# symmetric hump has ~zero linear correlation with date by construction,
# regardless of how strong the underlying relationship is. Test for it
# properly by adding a quadratic date term and asking whether it explains
# additional variance beyond the linear model, using the same trap-blocked
# permutation scheme as the linear date test (date permuted within trap,
# site/lure held fixed).

test_axis_date_quadratic <- function(axis_vals, dat_base, n_perm) {
  dat <- dat_base
  dat$y <- axis_vals

  fit_F <- function(d) {
    d$date_c <- as.numeric(d$date) - mean(as.numeric(d$date))  # center to avoid huge-magnitude collinear squared term
    m0 <- lm(y ~ site + lure + date_c, data = d)
    m1 <- lm(y ~ site + lure + date_c + I(date_c^2), data = d)
    anova(m0, m1)[2, "F"]
  }

  obs_F <- fit_F(dat)

  perm_F <- replicate(n_perm, {
    d2 <- dat %>% group_by(trap_id) %>% mutate(date = sample(date)) %>% ungroup()
    fit_F(d2)
  })

  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)
  c(F_stat = obs_F, p_perm = p_perm)
}

cat("\nTesting whether a quadratic date term explains additional variance beyond the linear model",
    "(", n_perm, "permutations each)...\n")

quad_results <- bind_rows(lapply(axis_names, function(ax) {
  r <- test_axis_date_quadratic(dat_base[[ax]], dat_base, n_perm)
  data.frame(axis = ax, F_stat = r["F_stat"], p_perm = r["p_perm"])
}))
rownames(quad_results) <- NULL
quad_results$q_value <- p.adjust(quad_results$p_perm, method = "BH")
cat("\n--- insect_PCoA1-3: quadratic date term, beyond linear (site+lure held fixed) ---\n")
print(quad_results, row.names = FALSE)
write.csv(quad_results, file.path(out_data_dir, "insect_pcoa_axes.quadratic_date_test.csv"), row.names = FALSE)

# Companion plot: loess fit (captures the hump shape) instead of a straight
# regression line, for visual comparison with Part 4's linear-fit version.
date_plot_loess <- ggplot(plot_dat, aes(x = date, y = value, color = site)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", color = "black", se = TRUE, linewidth = 0.5) +
  facet_wrap(~axis, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Collection date", y = "Axis value",
       title = "insect_PCoA1-3 vs. collection date (loess fit, prevalence-filtered insect table)") +
  theme_bw()
ggsave(file.path(out_fig_dir, "insect_pcoa_top3axes.by_date.loess.png"), date_plot_loess, width = 10, height = 4, dpi = 150)

## ---- 7. Which insect taxa drive each PCoA axis ----------------------------
# Pearson correlation of each insect taxon's Hellinger-transformed abundance
# (the exact data insect_PCoA1-3 were computed from -- Hellinger + Euclidean
# = Hellinger distance, see insect_fungal_coca_analysis.general_workflow.r)
# against each axis's scores. Same logic as the informal PCoA1-driver check
# from earlier in this project (e.g. "Melanophthalma sp., cor = 0.96" in
# iterative_analysis_updates.md), just formalized here and extended to all
# three axes and the full prevalence-filtered (all-families) taxon set.

insect_hel_df <- as.data.frame(insect_hel)
stopifnot(all(rownames(insect_hel_df) == dat_base$sample_id))

taxon_family <- sp_tab %>% distinct(Finest.ID, Family, Subfamily, Order)

axis_drivers <- bind_rows(lapply(axis_names, function(ax) {
  axis_vals <- dat_base[[ax]]
  cors <- sapply(insect_hel_df, function(col) cor(col, axis_vals))
  data.frame(axis = ax, taxon = names(cors), r = unname(cors)) %>%
    arrange(desc(abs(r)))
})) %>%
  left_join(taxon_family, by = c("taxon" = "Finest.ID"))

cat("\n--- Top 15 insect taxa by |correlation| with each PCoA axis (Hellinger-transformed abundance) ---\n")
for (ax in axis_names) {
  cat("\n", ax, ":\n", sep = "")
  print(axis_drivers %>% filter(axis == ax) %>% select(taxon, r, Family, Subfamily) %>% head(15),
        row.names = FALSE)
}
write.csv(axis_drivers, file.path(out_data_dir, "insect_pcoa_axes.taxon_drivers.csv"), row.names = FALSE)

# Focused figure: top 10 PCoA2 drivers, since that's this follow-up's
# specific question -- direction (sign of r) shown by fill color.
pcoa2_top_drivers <- axis_drivers %>% filter(axis == "insect_PCoA2") %>% slice_max(abs(r), n = 10)
pcoa2_driver_plot <- ggplot(pcoa2_top_drivers, aes(x = reorder(taxon, r), y = r, fill = r > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "#b2182b"), guide = "none") +
  labs(x = NULL, y = "Correlation with insect_PCoA2 (Hellinger-transformed abundance)",
       title = "Top 10 insect taxa driving insect_PCoA2") +
  theme_bw()
ggsave(file.path(out_fig_dir, "insect_pcoa2_top_taxon_drivers.png"), pcoa2_driver_plot, width = 7, height = 5, dpi = 150)

cat("\nDone. Outputs written to:\n",
    " ", file.path(out_data_dir, "insect_pcoa_axes.date_site_lure_association.csv"), "\n",
    " ", file.path(out_data_dir, "PCoA2_significant_taxa_taxonomy.csv"), "\n",
    " ", file.path(out_data_dir, "insect_pcoa_axes.quadratic_date_test.csv"), "\n",
    " ", file.path(out_data_dir, "insect_pcoa_axes.taxon_drivers.csv"), "\n",
    " ", file.path(out_fig_dir, "insect_pcoa_top3axes.by_date.png"), "\n",
    " ", file.path(out_fig_dir, "insect_pcoa_top3axes.by_date.loess.png"), "\n",
    " ", file.path(out_fig_dir, "insect_pcoa_top3axes.by_lure.png"), "\n",
    " ", file.path(out_fig_dir, "insect_pcoa2_top_taxon_drivers.png"), "\n")
