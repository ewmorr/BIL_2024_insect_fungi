##############################################################################
# Insect Community PCoA1-3 -- Association with Date, Site, and Lure
#
# insect_fungal_coca_analysis.general_workflow.r builds a Hellinger+Euclidean
# PCoA of the insect community and uses its first 3 axes (insect_PCoA1-3) as
# predictors of fungal abundance, noting in passing that insect_PCoA1 is
# "strongly collinear with collection date (r ~ 0.78)" and insect_PCoA3 is
# "strongly collinear with lure (partial R2 = 0.81, permutation q = 0.013)".
# The lure claim was formalized in insect_pcoa_lure_association.r (a
# broader screen across the first 10 axes); the date claim was never given
# its own script and only existed as narrative text in interpretation.md.
# This script consolidates both: it regresses exactly the 3 axes used
# elsewhere (insect_PCoA1-3) against date, site, and lure, so all three
# numbers can be looked up in one place and cited directly instead of
# reconstructed from a chat log.
#
# Two different permutation schemes are needed because date and lure vary
# at different levels of the design:
#   - Lure is a trap-level constant (one trap per site x lure combination),
#     so it's tested by permuting lure labels among the 3 traps WITHIN each
#     site (same scheme as insect_pcoa_lure_association.r).
#   - Date is fully crossed with trap (all 12 traps sampled on the same 6
#     collection dates), so it's tested by permuting date labels among the
#     6 samples WITHIN each trap -- the mirror image of the lure scheme,
#     and consistent with the within-trap permutation logic already used
#     for repeated-measures non-independence elsewhere in this project
#     (insect_fungal_coca_analysis.general_workflow.r Step 5).
# Site is reported for descriptive completeness (as in the original
# interpretation.md narrative) but only with a parametric partial F-test --
# with just 4 sites there's no sensible permutation unit for it (same
# caveat noted in insect_exploratory/insect_ords.R).
#
# Uses the same 69-sample matched insect/fungal dataset (Curculionidae +
# Latridiidae trap samples with paired ITS2 read data) as the CoCA script.
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

n_perm <- 999

# Partial R2 and semi-partial R2 for `term`, both derived from the same
# sequential ANOVA of the full model with `term` entered last, so ss_term is
# term's unique contribution controlling for the other two covariates
# (equivalent to R2_full - R2_reduced-without-term, i.e. what you'd get from
# regressing term on the other covariates, regressing y on those residuals,
# and squaring that correlation -- the two approaches coincide because
# ss_term here is already the marginal/Type II sum of squares).
#   - partial R2      = ss_term / (ss_term + ss_resid_full)
#                        "share of the variance NOT explained by the other
#                        terms that this term explains" -- denominator
#                        shrinks per term, so these do NOT sum to <=1 across
#                        terms (see interpretation.md discussion).
#   - semi-partial R2  = ss_term / ss_total
#                        "share of the TOTAL variance uniquely explained by
#                        this term" -- bounded by the full model's R2, so
#                        these sum to <=1 (exactly R2_full, if terms were
#                        orthogonal; less than that here since site/lure/date
#                        share some explanatory power).
term_r2s <- function(dat, term) {
  rhs <- c(setdiff(c("site", "lure", "date"), term), term)
  full_aov <- anova(lm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = dat))
  ss_term <- full_aov[term, "Sum Sq"]
  ss_resid <- full_aov["Residuals", "Sum Sq"]
  ss_total <- sum((dat$y - mean(dat$y))^2)
  c(r2_partial = ss_term / (ss_term + ss_resid), r2_semipartial = ss_term / ss_total)
}

# Date: permute date labels among the samples WITHIN each trap (mirror image
# of the lure test below), testing whether a trap's specific date-axis
# alignment carries information beyond site + lure.
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

# Lure: permute lure labels among the 3 traps WITHIN each site (same scheme
# as insect_pcoa_lure_association.r).
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

# Site: parametric partial F-test only -- descriptive, no permutation (see
# header note).
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

# BH correction within each term across the 3 axes (site's parametric p-values
# are descriptive only and not corrected, consistent with test_axis_site's
# "no permutation" caveat).
results <- results %>%
  group_by(term) %>%
  mutate(q_value = if (p_type[1] == "permutation") p.adjust(p_value, method = "BH") else NA_real_) %>%
  ungroup() %>%
  arrange(axis, term)

cat("\n--- insect_PCoA1-3 vs. date / site / lure ---\n")
print(as.data.frame(results), row.names = FALSE)
write.csv(results, "data/2024_insect_data/insect_pcoa_axes.date_site_lure_association.csv", row.names = FALSE)

## ---- 4. Diagnostic plots ---------------------------------------------------

plot_dat <- dat_base %>% select(all_of(axis_names), site, lure, date) %>%
  pivot_longer(all_of(axis_names), names_to = "axis", values_to = "value")

date_plot <- ggplot(plot_dat, aes(x = date, y = value, color = site)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  facet_wrap(~axis, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Collection date", y = "Axis value", title = "insect_PCoA1-3 vs. collection date") +
  theme_bw()
print(date_plot)
ggsave("figures/insect_pcoa_top3axes.by_date.png", date_plot, width = 10, height = 4, dpi = 150)

lure_plot <- ggplot(plot_dat, aes(x = lure, y = value, fill = lure)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, aes(shape = site), size = 2) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_wrap(~axis, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Lure", y = "Axis value", title = "insect_PCoA1-3 by lure") +
  theme_bw()
print(lure_plot)
ggsave("figures/insect_pcoa_top3axes.by_lure.png", lure_plot, width = 10, height = 4, dpi = 150)

cat("\nDone. Outputs written to:\n",
    "  data/2024_insect_data/insect_pcoa_axes.date_site_lure_association.csv\n",
    "  figures/insect_pcoa_top3axes.by_date.png\n",
    "  figures/insect_pcoa_top3axes.by_lure.png\n")
