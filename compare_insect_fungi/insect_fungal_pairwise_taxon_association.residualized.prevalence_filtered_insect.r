##############################################################################
# Follow-up to insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r:
# does the signal survive a stricter seasonal control on the insect-predictor
# side?  -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of insect_fungal_pairwise_taxon_association.residualized.r.
# Method is unchanged from that script; only the insect table construction is
# swapped to the Lineage-B one (all families, >=5-sample prevalence filter),
# matching insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r.
#
# Why this test (unchanged rationale from the parent): with LINEAR date in
# both the insect predictor's model and the fungal response's model, the
# fungal-side coefficient on insect_taxon is mathematically IDENTICAL whether
# insect_taxon enters raw or pre-residualized against site+date first
# (Frisch-Waugh-Lovell) -- so "residualize against linear date first" would
# just reproduce the main-grid result. For a genuinely stricter test the
# insect taxon needs to be stripped of MORE seasonal signal than a straight
# line can remove. There are only 6 distinct collection dates, so we use the
# strongest seasonal control the data supports without interpolating between
# rounds: date as a 6-level FACTOR when residualizing each insect taxon.
# Because linear date is in the span of the 6 date-dummy columns, this
# guarantees insect_taxon_resid is EXACTLY orthogonal to site and to date in
# ANY functional form representable by 6 levels -- so if the main-grid result
# was driven by residual (nonlinear / stepwise) seasonal structure the linear
# term missed, it should collapse here.
#
# Fungal-side model keeps site + lure + factor(date) + insect_taxon_resid:
# lure is not part of the residualization (matching the parent -- "date/site"
# only), so it stays a live covariate; site/date are retained too even though
# insect_taxon_resid is exactly orthogonal to them by construction (inert for
# its coefficient, but they still absorb real fungal variance and tighten the
# residual noise the predictor is tested against).
#
# Outputs go to the Lineage-B _prevalence_filtered_insect directory alongside
# the main-grid script's outputs, with a .date_site_residualized suffix.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect"
out_fig_dir <- "figures/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match (identical to the Lineage-B main-grid script) ------

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]

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
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))
cat("6 distinct collection dates:", length(unique(meta$date)) == 6, "\n")

## ---- 2. Filter and transform (identical to the Lineage-B main-grid script) -

keep_f <- colSums(fungal > 0) >= 5
fungal_f <- fungal[, keep_f, drop = FALSE]
keep_i <- colSums(insect > 0) >= 5
cat(sum(keep_f), "fungal taxa,", sum(keep_i),
    "insect taxa retained after prevalence filter (matched samples).\n")

insect_hel_full <- decostand(insect, method = "hellinger")
insect_hel <- insect_hel_full[, keep_i, drop = FALSE]

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)

stopifnot(all(rownames(fungal_clr) == meta$sample_id))
stopifnot(all(rownames(insect_hel) == meta$sample_id))

n_insect_taxa <- ncol(insect_hel)
n_fungal_taxa <- ncol(fungal_clr)

## ---- 3. Residualize each insect taxon against site + factor(date) ---------

meta$date_f <- factor(meta$date)
Z_resid <- model.matrix(~ site + date_f, data = meta)
cat("\nResidualizing", n_insect_taxa, "insect taxa against site + factor(date) (",
    ncol(Z_resid), "-column design, factor(date) has", nlevels(meta$date_f), "levels)...\n")

insect_mat <- as.matrix(insect_hel)
resid_fit <- lm.fit(Z_resid, insect_mat)
insect_resid <- resid_fit$residuals
colnames(insect_resid) <- colnames(insect_mat)
rownames(insect_resid) <- rownames(insect_mat)

chk_taxon <- colnames(insect_resid)[1]
cat("Sanity check (should be ~0):",
    "cor(resid, date_numeric) =", round(cor(insect_resid[, chk_taxon], as.numeric(meta$date)), 10),
    "\n")

## ---- 4. Vectorized per-pair permutation test, fungal-side model keeps -----
## site + lure + factor(date) + insect_taxon_resid

X_fixed <- model.matrix(~ site + lure + date_f, data = meta)
p_fixed <- ncol(X_fixed)
focal_idx <- p_fixed + 1
n <- nrow(X_fixed)

Y <- as.matrix(fungal_clr)

vectorized_t <- function(X, Y, focal_idx) {
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% crossprod(X, Y)
  resid <- Y - X %*% beta
  df <- nrow(X) - ncol(X)
  sigma2 <- colSums(resid^2) / df
  se_focal <- sqrt(sigma2 * XtX_inv[focal_idx, focal_idx])
  beta[focal_idx, ] / se_focal
}

n_perm <- 999
ctrl <- how(within = Within(type = "free"), blocks = meta$trap_id, nperm = n_perm)
perm_ids <- shuffleSet(n, control = ctrl)

t0 <- Sys.time()
pairwise_list <- vector("list", n_insect_taxa)
names(pairwise_list) <- colnames(insect_resid)

for (itax in colnames(insect_resid)) {
  x_focal <- insect_resid[, itax]
  X_obs <- cbind(X_fixed, insect_taxon = x_focal)

  chk <- tryCatch(solve(crossprod(X_obs)), error = function(e) NULL)
  if (is.null(chk)) {
    cat("SKIPPING", itax, "-- design matrix is singular.\n")
    next
  }

  obs_t <- vectorized_t(X_obs, Y, focal_idx)

  exceed <- integer(n_fungal_taxa)
  for (i in seq_len(n_perm)) {
    X_perm <- cbind(X_fixed, insect_taxon = x_focal[perm_ids[i, ]])
    t_perm <- vectorized_t(X_perm, Y, focal_idx)
    exceed <- exceed + (abs(t_perm) >= abs(obs_t))
  }
  p_perm <- (exceed + 1) / (n_perm + 1)

  pairwise_list[[itax]] <- data.frame(
    insect_taxon = itax,
    fungal_taxon = colnames(Y),
    t_stat = obs_t,
    p_perm = p_perm,
    row.names = NULL
  )
}
cat("Full grid permutation test (residualized) completed in",
    round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min\n")

full_grid_resid <- bind_rows(pairwise_list)

full_grid_resid <- full_grid_resid %>%
  group_by(insect_taxon) %>%
  mutate(q_value_within_insect = p.adjust(p_perm, method = "BH")) %>%
  ungroup() %>%
  mutate(q_value_global = p.adjust(p_perm, method = "BH"))

cat("\nTotal pairs tested:", nrow(full_grid_resid), "\n")
cat("Significant at q_value_within_insect < 0.10:", sum(full_grid_resid$q_value_within_insect < 0.10), "\n")
cat("Significant at q_value_global < 0.10:", sum(full_grid_resid$q_value_global < 0.10), "\n")

write.csv(full_grid_resid,
          file.path(out_data_dir, "fungal_insect_pairwise_full_grid.date_site_residualized.csv"),
          row.names = FALSE)

## ---- 5. Compare directly against the Lineage-B main grid -----------------

orig <- read.csv(file.path(out_data_dir, "fungal_insect_pairwise_full_grid.csv"))

per_insect_orig <- orig %>% group_by(insect_taxon) %>% summarise(n_sig_original = sum(q_value_within_insect < 0.10))
per_insect_resid <- full_grid_resid %>% group_by(insect_taxon) %>% summarise(n_sig_residualized = sum(q_value_within_insect < 0.10))

comparison <- full_join(per_insect_orig, per_insect_resid, by = "insect_taxon") %>%
  mutate(across(starts_with("n_sig"), ~ replace_na(., 0))) %>%
  arrange(desc(n_sig_original))

cat("\n--- Main grid vs. site/factor(date)-residualized: n significant fungal partners per insect taxon (top 20) ---\n")
print(head(as.data.frame(comparison), 20), row.names = FALSE)
write.csv(comparison, file.path(out_data_dir, "original_vs_residualized_comparison.csv"), row.names = FALSE)

# Focus on the insect taxa that carried the most fungal partners in the main
# grid -- data-driven here (unlike the Lineage-A parent, which hard-coded the
# 3 early-season Scolytinae it had already flagged), since the Lineage-B hit
# set isn't known in advance.
flagged <- comparison %>% filter(n_sig_original > 0) %>% slice_max(n_sig_original, n = 5) %>% pull(insect_taxon)
for (tx in flagged) {
  orig_hits <- orig %>% filter(insect_taxon == tx, q_value_within_insect < 0.10) %>% pull(fungal_taxon)
  resid_hits <- full_grid_resid %>% filter(insect_taxon == tx, q_value_within_insect < 0.10) %>% pull(fungal_taxon)
  overlap <- length(intersect(orig_hits, resid_hits))
  cat("\n", tx, ": main grid n=", length(orig_hits), " residualized n=", length(resid_hits),
      " overlap=", overlap,
      if (length(orig_hits) > 0) paste0(" (", round(100 * overlap / length(orig_hits), 1), "% of main-grid hits survive)") else "", "\n")
}

## ---- 6. Updated significant-hit table + plot -----------------------------

insect_taxonomy <- sp_tab %>%
  select(insect_taxon = Finest.ID, insect_order = Order, insect_family = Family,
         insect_subfamily = Subfamily, insect_genus = Genus, insect_species = Species) %>%
  distinct()
fungal_tax_df <- fungal_taxonomy %>%
  rownames_to_column("fungal_taxon") %>%
  select(fungal_taxon, fungal_kingdom = Kingdom, fungal_phylum = Phylum, fungal_class = Class,
         fungal_order = Order, fungal_family = Family, fungal_genus = Genus, fungal_species = Species)

sig_hits_resid <- full_grid_resid %>%
  filter(q_value_within_insect < 0.10) %>%
  left_join(insect_taxonomy, by = "insect_taxon") %>%
  left_join(fungal_tax_df, by = "fungal_taxon") %>%
  mutate(direction = ifelse(t_stat > 0, "positive", "negative")) %>%
  arrange(q_value_within_insect)
write.csv(sig_hits_resid,
          file.path(out_data_dir, "fungal_insect_pairwise_significant_hits.date_site_residualized.csv"),
          row.names = FALSE)

cat("\n", nrow(sig_hits_resid), "significant pairs after residualization (of", nrow(full_grid_resid), "tested)\n")
cat(length(unique(sig_hits_resid$fungal_taxon)), "unique fungal taxa,",
    length(unique(sig_hits_resid$insect_taxon)), "unique insect taxa involved\n")

plot_taxa <- comparison %>% filter(n_sig_original > 0 | n_sig_residualized > 0) %>% pull(insect_taxon)
if (length(plot_taxa) > 0) {
  comparison_plot <- comparison %>%
    pivot_longer(starts_with("n_sig"), names_to = "model", values_to = "n_sig") %>%
    mutate(model = recode(model, n_sig_original = "main grid (linear date)",
                           n_sig_residualized = "residualized (site+factor(date) removed)")) %>%
    filter(insect_taxon %in% plot_taxa)

  bar_plot <- ggplot(comparison_plot, aes(x = reorder(insect_taxon, n_sig), y = n_sig, fill = model)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(x = NULL, y = "Number of significantly associated fungal taxa (q<0.10)",
         title = "Main grid vs. site/factor(date)-residualized insect predictor (Lineage B)",
         fill = NULL) +
    theme_minimal()
  ggsave(file.path(out_fig_dir, "original_vs_residualized_comparison.png"),
         bar_plot, width = 8, height = max(4, 0.3 * length(plot_taxa)), dpi = 150, limitsize = FALSE)
} else {
  cat("\nNo insect taxa with significant partners in either model -- skipping comparison plot.\n")
}

cat("\nDone. Outputs written to", out_data_dir, ":\n",
    "  fungal_insect_pairwise_full_grid.date_site_residualized.csv\n",
    "  fungal_insect_pairwise_significant_hits.date_site_residualized.csv\n",
    "  original_vs_residualized_comparison.csv\n",
    "and", out_fig_dir, "/original_vs_residualized_comparison.png\n")
