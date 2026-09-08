##############################################################################
# Insect-Fungal Community Association Analysis
#
# Step 1: Co-correspondence analysis (CoCA) to identify which fungal taxa
#         load most strongly on axes shared with insect community composition.
# Step 2: Per-taxon permutation test of fungal abundance ~ insect community
#         axis, with FDR correction, run both without and with collection
#         date as a covariate, to see how much of the association survives
#         controlling for a shared seasonal trend. Also run without lure as
#         a sanity check (lure strongly structures the insect community, but
#         there's no biological reason to expect  fungi respond to it
#         independently, i.e., some taxa may be favored by ethanol, but all the  
#         lures contain ethanol) to confirm lure isn't partialing out real signal.
# Step 3: Direct per-taxon test of insect and fungal abundance ~ collection
#         date, to identify which specific taxa drive the seasonal signal
#         (interesting in its own right, independent of whether it implies
#         a direct insect-fungal relationship).
#
# Insect community here is restricted to Curculionidae + Latridiidae (the
# two dominant families in the 2024 NH trap catch, see insect_exploratory/
# insect_ords.R), matched to fungal ITS2 ASV read counts from the same trap/
# date via SequenceID. Samples are trapID x Lure combinations sampled
# repeatedly over collection dates, so samples are NOT independent within
# trap -- see Step 5 for how that's handled.
#
# Step 5 tests insect_PCoA1-3 as predictors (originally just PCoA1/PCoA2).
# insect_pcoa_lure_association.r screened all insect PCoA axes for
# association with lure and found PCoA3 (10.2% variance) far more strongly
# lure-associated (partial R2 = 0.81, permutation q = 0.013) than PCoA1 or
# PCoA2 -- a candidate explanation for why the fungal PERMANOVA detects a
# small lure effect (fungal_community_seasonality.r: 10.1%, p=0.05) that
# never showed up among the top CoCA-loading fungal taxa, which were only
# ever tested against PCoA1/PCoA2. PCoA4 was also nominally significant for
# lure but is left out here per follow-up scope.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("cocorresp", "vegan", "permute",
                    "dplyr", "tidyr", "tibble", "ggplot2")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

library(cocorresp)
library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(1)

## ---- 1. Load and match insect + fungal data ---------------------------------
# insect_species_tab.csv : full trap-catch table (see insect_exploratory/),
#   columns Class..Finest.ID + one column per sample (col_names format,
#   e.g. "C1.E..5.1").
# insect_community_metadata.csv : SequenceID, sampleID, col_names, trapID,
#   Site, Lure, CollectionDate -- links insect sample columns to the fungal
#   SequenceID used in ASV_tab.csv.
# ASV_tab.csv : fungal ITS2 ASV read counts, rows = SequenceID.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

# restrict to the two dominant families, as in insect_exploratory/insect_ords.R
sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID

# drop singleton taxa (present as a single individual total) and any
# resulting all-zero sample rows
sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]

# map insect sample columns (col_names) -> fungal SequenceID, and restrict to
# samples present in both the insect and fungal tables
id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)

shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full),
    "insect samples have matching fungal ASV data.\n")

id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

# Restrict to ASVs confirmed as Fungi -- ASV_tab.csv includes non-fungal
# (plant/animal/protist) and taxonomically unidentified ASVs; see
# ASVs_taxonomy.tsv (Kingdom column).
fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
cat(sum(is_fungus), "of", ncol(fungal), "ASVs confirmed Kingdom == k__Fungi (dropping",
    sum(!is_fungus), "non-fungal/unidentified ASVs).\n")
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

cat("Final matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (Curculionidae + Latridiidae),", ncol(fungal), "fungal ASVs.\n")

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))

## ---- 2. Filter and transform -----------------------------------------------

# Drop very rare fungal taxa before testing -- singletons/doubletons add
# multiple-testing burden without power to detect anything.
keep <- colSums(fungal > 0) >= 5              # present in >=5 of the samples
fungal_f <- fungal[, keep, drop = FALSE]
cat(sum(keep), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

# CoCA operates on Hellinger-transformed (or similarly scaled) community data
insect_hel <- decostand(insect, method = "hellinger")
fungal_hel <- decostand(fungal_f, method = "hellinger")

## ---- 3. Co-correspondence analysis -----------------------------------------
# Predictive CoCA: fungal composition is the response (y), insect composition
# is the predictor (x). (Symmetric CoCA -- method = "symmetric" -- is the
# alternative if you don't want to designate a predictor/response direction.)

coca_full <- coca(y = fungal_hel, x = insect_hel, method = "predictive")

# --- Decide how many PLS axes to retain ---
# cocorresp provides two complementary diagnostics for this: leave-one-out
# cross-validation (the package authors' preferred method, but slow) and a
# permutation test of each axis's significance. Inspect the printed output
# of both, then set n_axes_use by hand -- there isn't a single field to read
# this off of automatically.

max_axes_to_check <- min(3, ncol(insect_hel) - 1)

cat("\nLeave-one-out cross-validation (can take a while)...\n")
coca_cv <- crossval(y = fungal_hel, x = insect_hel, n.axes = max_axes_to_check)
print(summary(coca_cv))

cat("\nPermutation test of axis significance...\n")
coca_perm <- permutest(coca_full, permutations = 99, n.axes = max_axes_to_check)
print(coca_perm)

# Set based on the CV press statistics / permutation p-values printed above.
# Defaulting to 2, consistent with the cocorresp package's own worked examples.
n_axes_use <- 3

coca_mod <- coca(y = fungal_hel, x = insect_hel, method = "predictive", n.axes = n_axes_use)

# Fungal (response, Y-block) species scores on the retained CoCA axes = how
# strongly each fungal taxon loads on the insect-explained structure.
sc <- scores(coca_mod, display = "species", choices = seq_len(n_axes_use))
fungal_scores <- as.data.frame(sc$species$Y)
colnames(fungal_scores) <- paste0("CoCA", seq_len(n_axes_use))
fungal_scores$taxon <- rownames(fungal_scores)

# Rank by vector length across the retained axes (overall association strength,
# axis-1-only if n_axes_use == 1)
score_cols <- paste0("CoCA", seq_len(n_axes_use))
fungal_scores$coca_strength <- sqrt(rowSums(fungal_scores[, score_cols, drop = FALSE]^2))
fungal_scores <- fungal_scores %>% arrange(desc(coca_strength))

cat("\nTop 15 fungal taxa by CoCA loading strength:\n")
print(head(fungal_scores[, c("taxon", "coca_strength")], 15))

# Quick biplot for visual inspection (response block = fungal community)
pdf("figures/compare_insects_fungi_top3axes/coca_biplot.pdf", width = 7, height = 7)
plot(coca_mod, which = "response", type = "text",
     main = "Predictive CoCA: fungal (response) ~ insect (predictor)")
dev.off()

## ---- 4. Insect community axis for use as a predictor in Step 5 -------------
# Use an independent ordination of the insect data (NMDS or PCoA) rather than
# the CoCA axes themselves, so Step 5's test isn't circular with Step 3.

insect_dist <- vegdist(insect_hel, method = "euclidean")  # Hellinger + Euclidean ~ chord distance
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)
insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

## ---- 5. Per-taxon mixed-model / permutation test ---------------------------
# Test whether each fungal taxon's abundance tracks the insect PCoA axes,
# accounting for the repeated-measures structure (same trap sampled across
# dates) and for site/lure as blocking factors.
#
# NOTE ON RANDOM EFFECTS: with only 4 sites, a random intercept for site has
# too few levels for a stable variance estimate. Site and lure are therefore
# fit as FIXED blocking factors; non-independence from repeated sampling of
# the same trap is instead handled via a permutation test that only shuffles
# insect-axis values WITHIN trap (i.e. across collection dates for a given
# trap), preserving the trap-level structure under the null.
#
# NOTE ON DATE: insect_PCoA1 is strongly collinear with collection date
# (r = 0.78, F = 168, permutation p = 0.001 --
# compare_insect_fungi/insect_pcoa_date_lure_association.r) -- it is mostly
# a seasonal turnover axis (early-season Scolytinae-dominated catches vs.
# later-season Latridiidae-dominated catches). We run this test twice,
# WITHOUT and WITH date as a fixed covariate, and keep both results: the
# no-date version answers "does fungal abundance track insect community
# composition" (which may just reflect a shared seasonal trend); the
# date-adjusted version answers the stricter question of whether that
# association holds *beyond* a shared seasonal trend. Comparing the two is
# itself informative -- a taxon that drops out once date is added is (at
# least partly) riding the season, not necessarily tracking the insects
# directly.
#
# NOTE ON LURE / PCoA3: insect_PCoA3 is strongly collinear with lure
# (partial R2 = 0.81, permutation q = 0.013 --
# compare_insect_fungi/insect_pcoa_lure_association.r; restricted to just
# PCoA1-3, compare_insect_fungi/insect_pcoa_date_lure_association.r gets
# partial R2 = 0.81, q = 0.0045), the mirror image of PCoA1's relationship
# to date. Because "lure" is also in rhs_terms below,
# testing insect_PCoA3 WITH lure included asks a narrow question -- do
# fungi track the ~19% of PCoA3 that ISN'T lure -- while dropping lure
# (include_lure = FALSE) asks the broader question of whether fungi track
# PCoA3/lure at all. The gap between those two answers is the direct test
# of whether this axis is what's producing the small lure effect the
# fungal PERMANOVA detects (fungal_community_seasonality.r: 10.1%, p=0.05).
#
# NOTE ON TAXON COUNT: fungal_f has thousands of ASVs (prevalence filter
# alone isn't restrictive for a real ITS2 table), and each permutation test
# below is a full model refit x 999 permutations. Running that on every
# retained ASV is both computationally impractical (~1 hour+) and a poor use
# of the permutation test, which is meant to confirm the CoCA screen's
# top candidates, not stand in for it. We test only the top-ranked taxa by
# CoCA loading strength (Step 3).

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}

fungal_clr <- clr_transform(fungal_f)

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

n_candidates <- 200
candidate_taxa <- head(fungal_scores$taxon, n_candidates)
cat("\nRunning per-taxon permutation test on the top", length(candidate_taxa),
    "of", ncol(fungal_clr), "fungal taxa (by CoCA loading strength).\n")

n_perm <- 999

test_one_taxon <- function(taxon_abund, dat_base, n_perm, test_axis, include_date, include_lure = TRUE, block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  rhs_terms <- c("site", if (include_lure) "lure", if (include_date) "date", "insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
  model_formula <- as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))

  fit_stat <- function(d) {
    m <- lm(model_formula, data = d)
    cs <- coef(summary(m))
    c(t = cs[test_axis, "t value"], b = cs[test_axis, "Estimate"])
  }

  obs <- fit_stat(dat)
  obs_t <- unname(obs["t"])

  # Fully-standardized ("beta weight") partial slope of the axis under test,
  # from the OBSERVED fit only: b * sd(x) / sd(y) -- the QuantPsyc::lm.beta /
  # effectsize::standardize_parameters(method = "basic") convention. Used as
  # the fig3 volcano x-axis instead of the raw t-statistic. The permutation
  # null (p_perm / q_value) is UNCHANGED and still built from the
  # t-statistic; sd(x) is invariant under the within-trap axis permutation
  # (same values, reordered) and y is never permuted.
  beta_std <- unname(obs["b"]) * sd(dat[[test_axis]]) / sd(dat$y)

  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)

  # permute only the axis under test; hold the other PCoA axis (and all other
  # covariates) fixed, so the null is specific to that axis's association
  perm_t <- apply(perm_ids, 1, function(idx) {
    d2 <- dat
    d2[[test_axis]] <- dat[[test_axis]][idx]
    unname(fit_stat(d2)["t"])
  })

  p_perm <- (sum(abs(perm_t) >= abs(obs_t)) + 1) / (n_perm + 1)
  c(t_stat = obs_t, p_perm = p_perm, beta_std = beta_std)
}

run_pcoa_association_test <- function(test_axis, include_date, include_lure = TRUE) {
  res <- bind_rows(lapply(candidate_taxa, function(tax) {
    r <- test_one_taxon(fungal_clr[, tax], dat_base, n_perm, test_axis = test_axis,
                         include_date = include_date, include_lure = include_lure)
    data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"], beta_std = r["beta_std"])
  }))
  res$q_value <- p.adjust(res$p_perm, method = "BH")
  res %>%
    left_join(fungal_scores %>% select(taxon, coca_strength), by = "taxon") %>%
    arrange(q_value)
}

pcoa_test_specs <- list(
  list(axis = "insect_PCoA1", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_only.csv"),
  list(axis = "insect_PCoA1", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_plus_date.csv"),
  list(axis = "insect_PCoA2", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA2_only.csv"),
  list(axis = "insect_PCoA2", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA2_plus_date.csv"),
  list(axis = "insect_PCoA3", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA3_only.csv"),
  list(axis = "insect_PCoA3", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA3_plus_date.csv"),
  # Sanity check: lure is a strong driver of the INSECT community (see Step 4
  # below) but there is no biological reason to expect fungi respond to lure
  # independently (all the lures contain ethanol), and a companion dataset
  # showed minimal fungal~lure effect.
  # Drop lure from the PCoA1 model and compare -- if lure were partialing out
  # real insect-fungal signal, dropping it should recover more significant
  # hits; if it's just an inert covariate, results should barely change.
  list(axis = "insect_PCoA1", include_date = FALSE, include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_only.no_lure.csv"),
  list(axis = "insect_PCoA1", include_date = TRUE,  include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_plus_date.no_lure.csv"),
  # Same check for PCoA3, and more consequential here: insect_PCoA3 is
  # itself ~81% explained by lure (see NOTE ON LURE / PCoA3 above), far more
  # than PCoA1/PCoA2. With lure included, a PCoA3 hit means "fungi track
  # the ~19% of PCoA3 that isn't lure"; dropping lure tests whether fungi
  # track PCoA3/lure more broadly, including the part that IS lure.
  # Comparing the two is the direct test of whether insect_PCoA3 is what's
  # producing the small lure effect the fungal PERMANOVA finds.
  list(axis = "insect_PCoA3", include_date = FALSE, include_lure = FALSE, file = "fungal_insect_association_results.PCoA3_only.no_lure.csv"),
  list(axis = "insect_PCoA3", include_date = TRUE,  include_lure = FALSE, file = "fungal_insect_association_results.PCoA3_plus_date.no_lure.csv")
)

pcoa_results <- lapply(pcoa_test_specs, function(spec) {
  cat("\n---", spec$axis, "association", if (spec$include_date) "WITH" else "WITHOUT", "date covariate,",
      if (spec$include_lure) "WITH" else "WITHOUT", "lure covariate ---\n")
  res <- run_pcoa_association_test(test_axis = spec$axis, include_date = spec$include_date, include_lure = spec$include_lure)
  cat(sum(res$q_value < 0.10), "of", nrow(res), "taxa significant at q<0.10\n")
  print(head(res, 15))
  write.csv(res, file.path("data/compare_insects_fungi_top3axes", spec$file), row.names = FALSE)
  res
})
names(pcoa_results) <- sapply(pcoa_test_specs, `[[`, "file")

results_with_date <- pcoa_results[["fungal_insect_association_results.PCoA1_plus_date.csv"]]
results_pcoa3_with_date <- pcoa_results[["fungal_insect_association_results.PCoA3_plus_date.csv"]]

# Generalized so the same lure-partialling check can be run for both
# PCoA1 (the original sanity check -- lure structures the insects but
# shouldn't independently structure the fungi) and PCoA3 (the axis
# insect_pcoa_lure_association.r flagged as the one actually carrying the
# lure signal).
lure_partialling_check <- function(axis_label) {
  key <- function(suffix) sprintf("fungal_insect_association_results.%s%s.csv", axis_label, suffix)
  summary_df <- data.frame(
    model = paste0(axis_label, c(", with lure, no date", ", without lure, no date",
                                  ", with lure + date", ", without lure + date")),
    n_sig = c(sum(pcoa_results[[key("_only")]]$q_value < 0.10),
              sum(pcoa_results[[key("_only.no_lure")]]$q_value < 0.10),
              sum(pcoa_results[[key("_plus_date")]]$q_value < 0.10),
              sum(pcoa_results[[key("_plus_date.no_lure")]]$q_value < 0.10))
  )
  print(summary_df, row.names = FALSE)

  t_cor_nodate <- cor(pcoa_results[[key("_only")]]$t_stat,
                       pcoa_results[[key("_only.no_lure")]]$t_stat[
                         match(pcoa_results[[key("_only")]]$taxon,
                               pcoa_results[[key("_only.no_lure")]]$taxon)])
  cat("Correlation of", axis_label, "t-statistics (no-date models), with vs without lure:", round(t_cor_nodate, 4), "\n")
  invisible(summary_df)
}

cat("\n--- Lure partialling check: does including lure change the insect_PCoA1 signal? ---\n")
lure_partialling_check("PCoA1")

cat("\n--- Lure partialling check: does including lure change the insect_PCoA3 signal? ---\n")
lure_partialling_check("PCoA3")

# Direct check: do the candidate fungal taxa respond to lure at all, independent
# of the insect axes? lure is constant within trap (trap = site x lure), so
# within-trap permutation can't test it meaningfully here -- a parametric,
# site-controlled F-test is used instead as a quick descriptive check (not
# permutation-corrected, unlike the tests above).
cat("\n--- Direct fungal ~ lure check (parametric, site-controlled; descriptive only) ---\n")
lure_direct_check <- bind_rows(lapply(candidate_taxa, function(tax) {
  d <- dat_base
  d$y <- fungal_clr[, tax]
  m0 <- lm(y ~ site, data = d)
  m1 <- lm(y ~ site + lure, data = d)
  a <- anova(m0, m1)
  data.frame(taxon = tax, p_lure = a[["Pr(>F)"]][2])
}))
lure_direct_check$q_lure <- p.adjust(lure_direct_check$p_lure, method = "BH")
cat(sum(lure_direct_check$q_lure < 0.10, na.rm = TRUE), "of", nrow(lure_direct_check),
    "candidate fungal taxa show a nominally significant lure effect (q<0.10, parametric, not permutation-corrected)\n")
write.csv(lure_direct_check, "data/compare_insects_fungi_top3axes/fungal_taxa_lure_direct_check.csv", row.names = FALSE)

## ---- 6. Diagnostic plots (date-adjusted models, the stricter test) --------

make_volcano <- function(res, subtitle) {
  ggplot(res, aes(x = t_stat, y = -log10(q_value))) +
    geom_point(aes(color = q_value < 0.10)) +
    geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
    labs(x = paste0("t-statistic (fungal abundance ~ ", subtitle, ")"),
         y = "-log10(FDR q-value)",
         title = paste("Fungal taxa associated with", subtitle),
         color = "q < 0.10") +
    theme_minimal()
}

volcano <- make_volcano(results_with_date, "insect PCoA1, date-adjusted")
print(volcano)
ggsave("figures/compare_insects_fungi_top3axes/fungal_insect_volcano.png", volcano, width = 6, height = 5, dpi = 150)

# PCoA3, with vs without lure -- the direct visual companion to the lure
# partialling check above: does dropping lure from the model reveal a
# fungal signal on PCoA3 that was suppressed while lure was included?
volcano_pcoa3 <- make_volcano(results_pcoa3_with_date, "insect PCoA3, date-adjusted, with lure")
print(volcano_pcoa3)
ggsave("figures/compare_insects_fungi_top3axes/fungal_insect_volcano.PCoA3.with_lure.png", volcano_pcoa3, width = 6, height = 5, dpi = 150)

volcano_pcoa3_nolure <- make_volcano(pcoa_results[["fungal_insect_association_results.PCoA3_plus_date.no_lure.csv"]],
                                      "insect PCoA3, date-adjusted, without lure")
print(volcano_pcoa3_nolure)
ggsave("figures/compare_insects_fungi_top3axes/fungal_insect_volcano.PCoA3.no_lure.png", volcano_pcoa3_nolure, width = 6, height = 5, dpi = 150)

## ---- 7. Direct seasonal (date) association ----------------------------------
# Complementary to Step 5: rather than asking whether fungal abundance
# tracks insect community composition, ask the simpler question of which
# individual insect and fungal taxa vary with collection date at all. The
# seasonal trend is a real and interesting part of the story on its own
# (see insect_PCoA1's drivers), even where it doesn't imply a direct
# fungal-insect relationship.

test_date_taxon <- function(taxon_abund, dat_base, n_perm, block_var = "trap_id") {
  # Tests a linear date term (t_stat/p_perm/q_value -- unchanged from the
  # original version of this function; every downstream script that reads
  # this output's .csv filters on q_value, so these columns/semantics are
  # kept exactly as before) AND a quadratic term (date centered to avoid a
  # huge-magnitude collinear squared term) in the same model, so each taxon
  # can also be classified as a hump/dip (non-monotonic) seasonal pattern --
  # see the shape classification below. Both terms come from the same model
  # fit, so this adds negligible extra cost per permutation.
  dat <- dat_base
  dat$y <- taxon_abund

  fit_stats <- function(d) {
    d$date_c <- as.numeric(d$date) - mean(as.numeric(d$date))
    m <- lm(y ~ site + lure + date_c + I(date_c^2), data = d)
    cs <- coef(summary(m))
    c(t_linear = cs["date_c", "t value"], t_quad = cs["I(date_c^2)", "t value"],
      b_linear = cs["date_c", "Estimate"], b_quad = cs["I(date_c^2)", "Estimate"])
  }

  obs <- fit_stats(dat)

  # Fully-standardized ("beta weight") partial slopes from the OBSERVED fit
  # only: b * sd(x) / sd(y) -- the QuantPsyc::lm.beta /
  # effectsize::standardize_parameters(method = "basic") convention. Reported
  # as an effect-size companion to the t-statistic (used as the fig3 volcano
  # x-axis instead of the raw t-statistic); the permutation null below
  # (p_perm / q_value) is UNCHANGED and still built from the t-statistic.
  # sd(x) is invariant under the within-trap date permutation (same values,
  # reordered) and y is never permuted, so the observed-fit value is the
  # correct one. For the quadratic term x = date_c^2, so beta_std_quad is
  # "per SD of date_c^2": its SIGN is concavity (hump < 0 / dip > 0, same as
  # t_stat_quad), its magnitude is not on the same scale as beta_std.
  date_c_obs <- as.numeric(dat$date) - mean(as.numeric(dat$date))
  sd_y <- sd(dat$y)
  beta_std      <- unname(obs["b_linear"]) * sd(date_c_obs)   / sd_y
  beta_std_quad <- unname(obs["b_quad"])   * sd(date_c_obs^2) / sd_y

  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)

  perm_stats <- t(apply(perm_ids, 1, function(idx) {
    d2 <- dat
    d2$date <- dat$date[idx]
    fit_stats(d2)
  }))

  p_linear <- (sum(abs(perm_stats[, "t_linear"]) >= abs(obs["t_linear"])) + 1) / (n_perm + 1)
  p_quad <- (sum(abs(perm_stats[, "t_quad"]) >= abs(obs["t_quad"])) + 1) / (n_perm + 1)

  c(t_stat = unname(obs["t_linear"]), p_perm = p_linear,
    t_stat_quad = unname(obs["t_quad"]), p_perm_quad = p_quad,
    beta_std = beta_std, beta_std_quad = beta_std_quad)
}

classify_shape <- function(df) {
  df %>%
    mutate(shape = case_when(
      q_value_quad < 0.10 & t_stat_quad < 0 ~ "hump (peaks mid-season)",
      q_value_quad < 0.10 & t_stat_quad > 0 ~ "dip (troughs mid-season)",
      q_value < 0.10 & t_stat > 0 ~ "linear increase (late-season)",
      q_value < 0.10 & t_stat < 0 ~ "linear decrease (early-season)",
      TRUE ~ "no significant date pattern"
    )) %>%
    arrange(pmin(q_value, q_value_quad))
}

# --- Insect taxa: only ncol(insect_hel) taxa after filtering -- cheap to test all ---
insect_hel_df <- as.data.frame(insect_hel)
stopifnot(all(rownames(insect_hel_df) == meta$sample_id))

cat("\nTesting", ncol(insect_hel_df), "insect taxa for direct association with collection date (linear + quadratic)...\n")
insect_date_results <- bind_rows(lapply(colnames(insect_hel_df), function(tax) {
  r <- test_date_taxon(insect_hel_df[[tax]], meta, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"],
             t_stat_quad = r["t_stat_quad"], p_perm_quad = r["p_perm_quad"],
             beta_std = r["beta_std"], beta_std_quad = r["beta_std_quad"])
}))
insect_date_results$q_value <- p.adjust(insect_date_results$p_perm, method = "BH")
insect_date_results$q_value_quad <- p.adjust(insect_date_results$p_perm_quad, method = "BH")
insect_date_results <- classify_shape(insect_date_results)
cat(sum(insect_date_results$q_value < 0.10), "of", nrow(insect_date_results),
    "insect taxa significantly associated with date (q<0.10, LINEAR term)\n")
cat(sum(insect_date_results$q_value_quad < 0.10), "of", nrow(insect_date_results),
    "show a significant QUADRATIC date term (q<0.10) -- real hump/dip shape.\n")
cat("\nShape breakdown:\n")
print(table(insect_date_results$shape))
print(head(insect_date_results, 15))
write.csv(insect_date_results, "data/2024_insect_data/insect_taxa_date_association.csv", row.names = FALSE)

# --- Fungal taxa: screen by raw correlation with date first (cheap, all
# prevalence-filtered ASVs), then confirm the top candidates with the
# permutation test -- the same two-stage screen -> confirm logic as
# Step 3/5, just screening on date instead of CoCA loading. ---
date_numeric <- as.numeric(dat_base$date)
fungal_date_cor <- apply(fungal_clr, 2, function(x) cor(x, date_numeric))
fungal_date_candidates <- names(sort(abs(fungal_date_cor), decreasing = TRUE))[seq_len(n_candidates)]

cat("\nTesting top", length(fungal_date_candidates), "of", ncol(fungal_clr),
    "fungal taxa (by |correlation with date|) for direct association with collection date (linear + quadratic)...\n")
fungal_date_results <- bind_rows(lapply(fungal_date_candidates, function(tax) {
  r <- test_date_taxon(fungal_clr[, tax], dat_base, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"],
             t_stat_quad = r["t_stat_quad"], p_perm_quad = r["p_perm_quad"],
             beta_std = r["beta_std"], beta_std_quad = r["beta_std_quad"])
}))
fungal_date_results$q_value <- p.adjust(fungal_date_results$p_perm, method = "BH")
fungal_date_results$q_value_quad <- p.adjust(fungal_date_results$p_perm_quad, method = "BH")
fungal_date_results <- classify_shape(fungal_date_results)
cat(sum(fungal_date_results$q_value < 0.10), "of", nrow(fungal_date_results),
    "tested fungal taxa significantly associated with date (q<0.10, LINEAR term)\n")
cat(sum(fungal_date_results$q_value_quad < 0.10), "of", nrow(fungal_date_results),
    "show a significant QUADRATIC date term (q<0.10) -- real hump/dip shape.\n")
cat("\nShape breakdown:\n")
print(table(fungal_date_results$shape))
print(head(fungal_date_results, 15))
write.csv(fungal_date_results, "data/compare_insects_fungi_top3axes/fungal_taxa_date_association.csv", row.names = FALSE)

cat("\nDone. Outputs written to data/compare_insects_fungi_top3axes/:\n",
    "  fungal_insect_association_results.PCoA1_only.csv\n",
    "  fungal_insect_association_results.PCoA1_plus_date.csv\n",
    "  fungal_insect_association_results.PCoA2_only.csv\n",
    "  fungal_insect_association_results.PCoA2_plus_date.csv\n",
    "  fungal_insect_association_results.PCoA3_only.csv\n",
    "  fungal_insect_association_results.PCoA3_plus_date.csv\n",
    "  fungal_insect_association_results.PCoA1_only.no_lure.csv\n",
    "  fungal_insect_association_results.PCoA1_plus_date.no_lure.csv\n",
    "  fungal_insect_association_results.PCoA3_only.no_lure.csv\n",
    "  fungal_insect_association_results.PCoA3_plus_date.no_lure.csv\n",
    "  fungal_taxa_lure_direct_check.csv\n",
    "  fungal_taxa_date_association.csv\n",
    "and data/2024_insect_data/:\n",
    "  insect_taxa_date_association.csv\n",
    "and figures/compare_insects_fungi_top3axes/:\n",
    "  coca_biplot.pdf\n",
    "  fungal_insect_volcano.png (PCoA1)\n",
    "  fungal_insect_volcano.PCoA3.with_lure.png\n",
    "  fungal_insect_volcano.PCoA3.no_lure.png\n")
