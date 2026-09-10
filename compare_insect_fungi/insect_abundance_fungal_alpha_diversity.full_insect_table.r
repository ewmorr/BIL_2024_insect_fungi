##############################################################################
# Total insect ABUNDANCE vs. fungal alpha diversity -- Lineage C (full insect
# table), new-hypothesis add-on to insect_fungal_alpha_diversity.full_insect_
# table.r.
#
# Hypothesis being tested
# -----------------------
# The main Lineage C script asks whether insect and fungal alpha DIVERSITY
# covary (i.e. whether diverse insect assemblages carry diverse fungal
# assemblages -- a pattern you'd expect if insect-fungus associations are
# partner-specific). It finds that link is weak (only Shannon, rho=-0.27,
# p=0.03; richness and Simpson null -- see the interp doc section 2/5).
#
# An alternative mechanism: insects act as NON-TARGETED dispersal vectors,
# picking up whatever fungal propagules they contact as they move through the
# environment. Under that model the relevant insect quantity is not how many
# DIFFERENT insect taxa a trap caught but how many INDIVIDUALS it caught --
# more insect traffic = more independent draws from the environmental fungal
# pool = more fungal taxa detected. Prediction: total insect abundance (sum
# of individuals across all taxa, per sample) correlates POSITIVELY with
# fungal alpha diversity, most clearly with fungal richness.
#
# Design
# ------
# - Insect abundance: rowSums() of the SAME raw all-families 426-taxon
#   trap-catch table the main Lineage C script uses (raw counts, empty-sample
#   drop only, NO singleton/prevalence filter -- project_organization.md's
#   "When to prevalence-filter vs. not"). One number per sample:
#   total_individuals. The per-sample insect richness is recomputed here too
#   purely as a cross-check that this table matches the main script's.
# - Fungal alpha diversity: NOT recomputed -- read straight from
#   fungal_alpha_diversity.csv (Shannon / Simpson dominance / richness, each
#   rarefied to 5000 x 100 iterations and averaged) written by
#   insect_fungal_alpha_diversity.full_insect_table.r, so the two scripts
#   can't drift. The asymptotic (Chao) fungal estimates from
#   insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r are also
#   read, if present, and the same correlations reported against them as a
#   secondary block (matching the observed-plus-asymptotic parallel the
#   interp doc keeps throughout).
#
# Two analyses per fungal metric (the seasonality caveat is the whole point):
#   1. raw_pooled -- Spearman(total_individuals, fungal_metric) over all
#      matched samples. Same plain unconstrained Spearman as the main
#      script's section 5 cross-community correlation. BUT both total insect
#      abundance and every fungal alpha metric have strong collection-date
#      structure in this dataset (interp doc section 4), so a raw pooled
#      correlation cannot by itself distinguish "insect traffic disperses
#      fungi" from "both just peak at the same point in the season".
#   2. season_residualized -- residualize log1p(total_individuals) AND the
#      fungal metric on `site + factor(date)` (the exact deseasonalising
#      convention of this project's *.residualized.r scripts: site + a
#      6-level date factor that absorbs ANY functional form of the seasonal
#      trend), then Spearman-correlate the residuals. This is the test that
#      actually bears on the dispersal-vector hypothesis: does a trap that
#      caught more insects THAN EXPECTED FOR ITS SITE AND DATE also yield
#      more fungal diversity than expected for its site and date? p-value by
#      free permutation of one residual vector (residuals are exchangeable
#      once site/date are removed); the parametric Spearman p on residuals is
#      also reported but is mildly anticonservative (ignores the df spent on
#      residualising) so the permutation p is the one to read.
#
# 3 metrics x 2 analyses = 6 tests (x2 again for the asymptotic block). Small
# fixed a priori set, reported with uncorrected p-values -- same convention
# as the site/lure/date covariate tables elsewhere in this project, not a
# taxon screen needing BH-FDR.
##############################################################################

## ---- 0. Setup -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_alpha_diversity_full_insect_table"
out_fig_dir  <- "figures/compare_insects_fungi_alpha_diversity_full_insect_table"
asym_data_dir <- "data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"
stopifnot(dir.exists(out_data_dir), dir.exists(out_fig_dir))

n_perm <- 9999

## ---- 1. Total insect abundance from the raw all-families table ---------------
# Reconstructed exactly as in insect_fungal_alpha_diversity.full_insect_table.r
# sections 1-3 (raw counts, drop empty samples only), then matched to the
# fungal alpha-diversity samples.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID
insect_full <- sp_tab.t[rowSums(sp_tab.t) > 0, ]                 # empty-sample drop only
cat("Raw insect table:", nrow(insect_full), "samples x", ncol(insect_full), "taxa (no filter).\n")

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))
insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID

insect_abund <- data.frame(
  sample_id        = rownames(insect),
  total_individuals = rowSums(insect),
  insect_richness_check = vegan::specnumber(insect),
  row.names = NULL
)

## ---- 2. Read fungal alpha diversity (do NOT recompute) ----------------------

fungal_alpha <- read.csv(file.path(out_data_dir, "fungal_alpha_diversity.csv"))
insect_alpha <- read.csv(file.path(out_data_dir, "insect_alpha_diversity.csv"))

# Cross-check: the insect table reconstructed here must match the one the main
# script used (same richness per sample).
chk <- insect_abund %>%
  inner_join(insect_alpha %>% select(sample_id, richness_main = richness), by = "sample_id")
if (!isTRUE(all.equal(chk$insect_richness_check, chk$richness_main))) {
  stop("Reconstructed insect richness does not match insect_alpha_diversity.csv -- ",
       "insect table construction has drifted from the main Lineage C script.")
}
cat("Insect-table cross-check OK (per-sample richness matches the main script).\n")

dat <- fungal_alpha %>%
  rename(fungal_shannon = shannon, fungal_simpson_dominance = simpson_dominance,
         fungal_richness = richness) %>%
  inner_join(insect_abund %>% select(sample_id, total_individuals), by = "sample_id") %>%
  inner_join(insect_alpha %>% select(sample_id, insect_shannon = shannon,
                                     insect_simpson_dominance = simpson_dominance,
                                     insect_richness = richness),
             by = "sample_id") %>%
  mutate(date = as.Date(date), log_total_individuals = log1p(total_individuals))
cat("\nMatched:", nrow(dat), "samples with total insect abundance + fungal alpha diversity.\n")
cat("Total individuals per sample: min", min(dat$total_individuals),
    "median", median(dat$total_individuals), "max", max(dat$total_individuals), "\n")

write.csv(dat %>% select(sample_id, site, lure, date, trap_id, total_individuals,
                         insect_shannon, insect_simpson_dominance, insect_richness,
                         fungal_shannon, fungal_simpson_dominance, fungal_richness),
          file.path(out_data_dir, "insect_abundance_by_sample.csv"), row.names = FALSE)

## ---- 3. Correlation engine -------------------------------------------------

fungal_metrics <- c("shannon", "simpson_dominance", "richness")
metric_labels  <- c(shannon = "Fungal Shannon diversity",
                    simpson_dominance = "Fungal Simpson dominance",
                    richness = "Fungal richness (rarefied)")

# raw pooled Spearman + free-permutation p
raw_corr <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  obs <- unname(ct$estimate)
  perm <- replicate(n_perm, suppressWarnings(
    cor(sample(x), y, method = "spearman")))
  data.frame(rho = obs, p_spearman = ct$p.value,
             p_perm = (sum(abs(perm) >= abs(obs)) + 1) / (n_perm + 1))
}

# residualize both sides on site + factor(date), Spearman on residuals,
# free-permutation p (residuals exchangeable once site/date removed)
resid_corr <- function(xvar, yvar, d) {
  rx <- residuals(lm(reformulate(c("site", "factor(date)"), xvar), data = d))
  ry <- residuals(lm(reformulate(c("site", "factor(date)"), yvar), data = d))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman"))
  obs <- unname(ct$estimate)
  perm <- replicate(n_perm, suppressWarnings(
    cor(sample(rx), ry, method = "spearman")))
  data.frame(rho = obs, p_spearman = ct$p.value,
             p_perm = (sum(abs(perm) >= abs(obs)) + 1) / (n_perm + 1))
}

run_block <- function(d, fungal_prefix, block_label) {
  bind_rows(lapply(fungal_metrics, function(m) {
    ycol <- paste0(fungal_prefix, m)
    r1 <- raw_corr(d$total_individuals, d[[ycol]])
    r2 <- resid_corr("log_total_individuals", ycol, d)
    bind_rows(
      data.frame(fungal_value = block_label, metric = m, analysis = "raw_pooled",
                 predictor = "total_individuals", r1, n = nrow(d)),
      data.frame(fungal_value = block_label, metric = m, analysis = "season_residualized",
                 predictor = "log1p(total_individuals) resid ~ site+factor(date)", r2, n = nrow(d))
    )
  }))
}

## ---- 4. Observed (rarefied) fungal alpha diversity ------------------------

results_obs <- run_block(dat, "fungal_", "observed_rarefied")

# Context: how seasonal / richness-linked is total insect abundance itself?
context <- bind_rows(
  data.frame(comparison = "total_individuals ~ date (numeric)",
             rho = suppressWarnings(cor(dat$total_individuals, as.numeric(dat$date), method = "spearman"))),
  data.frame(comparison = "total_individuals ~ insect_richness",
             rho = suppressWarnings(cor(dat$total_individuals, dat$insect_richness, method = "spearman"))),
  data.frame(comparison = "total_individuals ~ insect_shannon",
             rho = suppressWarnings(cor(dat$total_individuals, dat$insect_shannon, method = "spearman")))
)
cat("\n--- Context: total insect abundance vs. date and vs. insect alpha diversity (Spearman rho) ---\n")
print(context, row.names = FALSE)

## ---- 5. Asymptotic (Chao) fungal diversity, if available -----------------

results_asym <- NULL
asym_file <- file.path(asym_data_dir, "fungal_asymptotic_diversity.csv")
if (file.exists(asym_file)) {
  fungal_asym <- read.csv(asym_file)
  # column names in that file: shannon/simpson/richness are Hill numbers
  # q=1/q=2/q=0 ("simpson" there is a DIVERSITY, higher = more even -- the
  # OPPOSITE direction from the observed "simpson_dominance" above; kept as-is
  # and labelled so downstream, not sign-flipped).
  acols <- names(fungal_asym)
  amap <- c(shannon = grep("^shannon", acols, value = TRUE)[1],
            simpson_dominance = grep("^simpson", acols, value = TRUE)[1],
            richness = grep("^richness", acols, value = TRUE)[1])
  fungal_asym2 <- fungal_asym %>%
    transmute(sample_id,
              fungalasym_shannon = .data[[amap["shannon"]]],
              fungalasym_simpson_dominance = .data[[amap["simpson_dominance"]]],
              fungalasym_richness = .data[[amap["richness"]]])
  dat_asym <- dat %>% select(sample_id, site, date, total_individuals, log_total_individuals) %>%
    inner_join(fungal_asym2, by = "sample_id")
  cat("\nAsymptotic block:", nrow(dat_asym), "samples.\n")
  results_asym <- run_block(dat_asym, "fungalasym_", "asymptotic_chao")
  results_asym$note <- "asymptotic 'simpson_dominance' col is Hill q=2 DIVERSITY (higher=more even), opposite sign to observed"
}

## ---- 6. Assemble + write ------------------------------------------------------

results <- bind_rows(results_obs, results_asym)
results <- results %>%
  mutate(across(c(rho, p_spearman, p_perm), ~round(.x, 4)))
cat("\n=== Total insect abundance vs. fungal alpha diversity ===\n")
print(as.data.frame(results), row.names = FALSE)

write.csv(results, file.path(out_data_dir, "insect_abundance_fungal_alpha_correlation.csv"),
          row.names = FALSE)
write.csv(context, file.path(out_data_dir, "insect_abundance_context_correlations.csv"),
          row.names = FALSE)

## ---- 7. Figure: 2 rows (raw / residualized) x 3 fungal metrics -------------

theme_ab <- theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        legend.position = "bottom")

panel_raw <- function(m) {
  ycol <- paste0("fungal_", m)
  rr <- results_obs %>% filter(metric == m, analysis == "raw_pooled")
  lab <- sprintf("rho = %.2f\np_perm = %.3f", rr$rho, rr$p_perm)
  ggplot(dat, aes(total_individuals, .data[[ycol]], fill = as.numeric(date))) +
    geom_point(shape = 21, size = 2.4, color = "black", stroke = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.5) +
    scale_x_log10() +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         midpoint = median(as.numeric(dat$date)), name = "Date") +
    annotate("label", x = min(dat$total_individuals), y = Inf, hjust = 0, vjust = 1.2,
             label = lab, size = 3, fill = alpha("white", 0.75)) +
    labs(x = "Total insect individuals (log scale)", y = metric_labels[m],
         subtitle = "raw pooled") +
    theme_ab
}

panel_resid <- function(m) {
  ycol <- paste0("fungal_", m)
  rx <- residuals(lm(reformulate(c("site", "factor(date)"), "log_total_individuals"), data = dat))
  ry <- residuals(lm(reformulate(c("site", "factor(date)"), ycol), data = dat))
  pd <- data.frame(rx, ry, date = dat$date)
  rr <- results_obs %>% filter(metric == m, analysis == "season_residualized")
  lab <- sprintf("rho = %.2f\np_perm = %.3f", rr$rho, rr$p_perm)
  ggplot(pd, aes(rx, ry, fill = as.numeric(date))) +
    geom_point(shape = 21, size = 2.4, color = "black", stroke = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.5) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         midpoint = median(as.numeric(dat$date)), name = "Date") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, label = lab,
             size = 3, fill = alpha("white", 0.75)) +
    labs(x = "Insect abundance residual (~ site + date)",
         y = paste0(metric_labels[m], " residual"),
         subtitle = "season-residualized") +
    theme_ab
}

panels <- c(lapply(fungal_metrics, panel_raw), lapply(fungal_metrics, panel_resid))
fig <- gridExtra::arrangeGrob(grobs = panels, nrow = 2,
         top = "Total insect abundance vs. fungal alpha diversity (Lineage C, full insect table)")
ggsave(file.path(out_fig_dir, "insect_abundance_vs_fungal_alpha.png"), fig,
       width = 13, height = 8, dpi = 150)
ggsave(file.path(out_fig_dir, "insect_abundance_vs_fungal_alpha.pdf"), fig,
       width = 13, height = 8)

cat("\nDone. Outputs:\n",
    file.path(out_data_dir, "insect_abundance_fungal_alpha_correlation.csv"), "\n ",
    file.path(out_data_dir, "insect_abundance_by_sample.csv"), "\n ",
    file.path(out_data_dir, "insect_abundance_context_correlations.csv"), "\n ",
    file.path(out_fig_dir, "insect_abundance_vs_fungal_alpha.{png,pdf}"), "\n")
