##############################################################################
# Insect Community PCoA -- Which Axes Track Lure?
#
# insect_ords.R found roughly equal Site and Lure effects on the insect
# community by PERMANOVA (whole Bray-Curtis matrix). The CoCA workflow
# (insect_fungal_coca_analysis.general_workflow.r) then built a
# Hellinger+Euclidean PCoA of the same insect community and tested fungal
# taxa against only its first two axes as predictors -- PCoA1 turned out to
# be almost entirely a date/season axis (r=0.78 with collection date), and
# PCoA2 showed no association with anything. Lure explains a real (if
# smaller) share of the fungal PERMANOVA (10.1%, p=0.05,
# fungal_community_seasonality.r) despite no lure-associated taxa turning up
# in the CoCA loadings -- one explanation is that the lure signal in the
# insect community lives on a lower-variance axis (3rd, 4th, ...) that
# PCoA1/PCoA2 alone wouldn't capture, and that axis was never tested against
# fungal composition.
#
# This script does not touch the fungal data. It only asks: across the full
# set of insect PCoA axes (not just the first two), which ones are actually
# associated with lure? That determines whether it's worth extending the
# CoCA per-taxon tests (Step 5 of the CoCA script) to additional axes.
#
# Uses the same 69-sample matched insect/fungal dataset (Curculionidae +
# Latridiidae trap samples with paired ITS2 read data) as the CoCA script,
# so any axis flagged here can be plugged directly into that script's
# per-taxon test as an additional predictor.
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

## ---- 1. Load and match insect + fungal data (same as the CoCA script) -----
# Fungal data is loaded only to restrict to the same sample set used in the
# CoCA script -- fungal read counts themselves are never used below.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID

sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]

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
    "insect taxa (Curculionidae + Latridiidae).\n")

# Design check: this script's lure test relies on lure being constant within
# trap (one trap per site x lure combination) -- confirm that holds before
# using it as the permutation unit below.
trap_lure_map <- meta %>% distinct(trap_id, site, lure)
stopifnot(all(table(trap_lure_map$site, trap_lure_map$lure) == 1))
cat(nrow(trap_lure_map), "traps, one per site x lure combination, as expected.\n")

## ---- 2. Full insect PCoA (same Hellinger + Euclidean strategy as the CoCA
## script), retaining all non-trivial axes instead of just the top 2 --------
# Hellinger transform + Euclidean distance = Hellinger distance, one of the
# transformations Legendre & Gallagher (2001) recommend as a Bray-Curtis
# substitute for methods needing true Euclidean geometry (PCA, RDA, classical
# PCoA via cmdscale) -- see insect_fungal_coca_analysis.general_workflow.r
# for the same reasoning applied to just PCoA1/PCoA2.

insect_hel <- decostand(insect, method = "hellinger")
insect_dist <- vegdist(insect_hel, method = "euclidean")

max_k <- nrow(insect) - 1
insect_pcoa <- cmdscale(insect_dist, k = max_k, eig = TRUE)

cat("\nNegative eigenvalues:", sum(insect_pcoa$eig < 0),
    "of", length(insect_pcoa$eig),
    "(small/near-zero is expected numerical noise for Hellinger+Euclidean; large negative values would signal a non-Euclidean embedding)\n")

pos_eig <- insect_pcoa$eig[insect_pcoa$eig > 0]
var_explained <- pos_eig / sum(pos_eig)
cat("Variance explained by first 10 positive-eigenvalue axes:\n")
print(round(head(var_explained, 10), 4))

insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- paste0("insect_PCoA", seq_len(ncol(insect_axes)))
insect_axes$sample_id <- rownames(insect)

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(insect)))

## ---- 3. Per-axis lure association test --------------------------------------
# Lure is a trap-level constant (exactly one trap per site x lure
# combination, sampled repeatedly over dates -- see design check above), so
# it cannot be tested with the within-trap permutation used elsewhere in
# this project for taxon-level associations (that scheme permutes dates
# within a trap, which never changes a trap's lure). Instead, permute lure
# labels among the 3 traps WITHIN each site -- this preserves the site
# structure and the one-trap-per-lure design, and is the permutation that
# actually matches how lure varies in this experiment (across traps within
# a site, not across samples within a trap). With 3 lures x 4 sites, this
# is 6 possible within-site relabelings per site (3!), 6^4 = 1296 total
# distinct joint relabelings -- n_perm draws below are with replacement
# from that space.

n_axes_test <- min(10, ncol(insect_axes) - 1)
n_perm <- 999

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

  # descriptive partial R2 for lure (site, date already accounted for)
  ss <- anova(lm(y ~ site + date, data = dat), lm(y ~ site + lure + date, data = dat))
  r2_lure <- NA
  full_aov <- anova(lm(y ~ site + lure + date, data = dat))
  ss_resid <- full_aov["Residuals", "Sum Sq"]
  ss_lure <- full_aov["lure", "Sum Sq"]
  r2_lure <- ss_lure / (ss_lure + ss_resid)

  c(F_stat = obs_F, r2_lure = r2_lure, p_perm = p_perm)
}

cat("\nTesting lure association for the top", n_axes_test, "insect PCoA axes",
    "(", n_perm, "trap-within-site permutations each)...\n")

axis_names <- paste0("insect_PCoA", seq_len(n_axes_test))
lure_axis_results <- bind_rows(lapply(axis_names, function(ax) {
  r <- test_axis_lure(dat_base[[ax]], dat_base, trap_lure_map, n_perm)
  data.frame(axis = ax, pct_var_explained = 100 * var_explained[as.integer(sub("insect_PCoA", "", ax))],
             F_stat = r["F_stat"], r2_lure = r["r2_lure"], p_perm = r["p_perm"])
}))
lure_axis_results$q_value <- p.adjust(lure_axis_results$p_perm, method = "BH")
lure_axis_results <- lure_axis_results %>% arrange(p_perm)

cat("\n--- Insect PCoA axes ranked by lure association ---\n")
print(lure_axis_results, row.names = FALSE)
write.csv(lure_axis_results, "data/2024_insect_data/insect_pcoa_axes.lure_association.csv", row.names = FALSE)

top_lure_axis <- lure_axis_results$axis[1]
cat("\nStrongest lure-associated axis:", top_lure_axis,
    "(p_perm =", lure_axis_results$p_perm[1], ", q =", round(lure_axis_results$q_value[1], 4),
    ", partial R2 =", round(lure_axis_results$r2_lure[1], 3), ")\n")

## ---- 4. Diagnostic plots ----------------------------------------------------

var_df <- data.frame(axis = paste0("PCoA", seq_along(var_explained)),
                      axis_num = seq_along(var_explained),
                      pct_var = 100 * var_explained) %>%
  filter(axis_num <= n_axes_test) %>%
  left_join(lure_axis_results %>% transmute(axis_num = as.integer(sub("insect_PCoA", "", axis)), q_value),
            by = "axis_num")

scree_plot <- ggplot(var_df, aes(x = reorder(axis, axis_num), y = pct_var, fill = q_value < 0.10)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey70"), na.value = "grey70",
                     name = "lure q<0.10") +
  labs(x = "Insect PCoA axis", y = "% variance explained",
       title = "Insect PCoA axes: variance explained vs. lure association") +
  theme_bw()
print(scree_plot)
ggsave("figures/insect_pcoa_scree.lure_highlighted.png", scree_plot, width = 7, height = 5, dpi = 150)

top_axis_plot_dat <- dat_base %>% select(all_of(top_lure_axis), lure, site, date)
colnames(top_axis_plot_dat)[1] <- "axis_value"

top_axis_plot <- ggplot(top_axis_plot_dat, aes(x = lure, y = axis_value, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, aes(shape = site), size = 2) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Lure", y = top_lure_axis,
       title = paste0(top_lure_axis, " by lure (strongest lure-associated axis)")) +
  theme_bw()
print(top_axis_plot)
ggsave("figures/insect_pcoa_top_lure_axis.by_lure.png", top_axis_plot, width = 7, height = 5, dpi = 150)

cat("\nDone. Outputs written to:\n",
    "  data/2024_insect_data/insect_pcoa_axes.lure_association.csv\n",
    "  figures/insect_pcoa_scree.lure_highlighted.png\n",
    "  figures/insect_pcoa_top_lure_axis.by_lure.png\n")
