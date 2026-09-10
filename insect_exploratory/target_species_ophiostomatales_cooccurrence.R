##############################################################################
# Presence/absence co-occurrence:
#   target bark-beetle species  x  Ophiostomatales fungi
#
# The three species the EDRR traps target -- Ips grandicollis, Ips pini,
# Dendroctonus valens (see target_genus_phenology.R) -- are bark beetles whose
# natural history is tied to ophiostomatoid ("blue-stain") fungi. This script
# asks, purely on shared presence/absence across the 69 samples that have both
# insect and fungal data, whether any Ophiostomatales fungal taxon co-occurs
# with any of the three beetle species more (or less) than expected by chance.
#
# Comparisons are beetle x fungal-taxon ONLY. No fungus-vs-fungus and no
# beetle-vs-beetle pairs are computed.
#
# Metric: Jaccard similarity   J = a / (a + b + c)
#   a = samples where BOTH present, b = beetle only, c = fungus only.
#   Shared absence (d) is excluded on purpose -- co-presence is the focus and
#   joint absence is uninformative here (most samples lack most rare taxa).
#
# Fungal taxon is rolled up three ways, each run as its own grid; a
# genus/species/order unit is "present" in a sample if ANY constituent
# Ophiostomatales ASV is present:
#   (1) species  -- named Ophiostomatales species, prevalence >= MIN_PREV  [headline]
#   (2) genus    -- the 8 Ophiostomatales genera + "Any Ophiostomatales"
#   (3) ASV      -- individual Ophiostomatales ASVs, prevalence >= MIN_PREV
#
# Significance: restricted permutation null. The beetle presence/absence
# vector is shuffled N_PERM times WITHIN trap (permute::how with blocks = site,
# plots = trap type "none", within = free -- the same site/trap structure the
# rest of this project's permutation tests preserve; identical call to
# insect_ords.prevalence_filter.R). Jaccard is recomputed for every pair each
# permutation.
#   p_pos = add-one Pr(J_perm >= J_obs)   (co-occurrence)
#   p_neg = add-one Pr(J_perm <= J_obs)   (avoidance)
#   p_perm = 2 * min(p_pos, p_neg), capped at 1; BH-FDR (q_perm) within grid.
# Fisher's exact test on the 2x2 is reported alongside as a structure-blind
# cross-check (q_fisher = BH within grid).
# z_score = (J_obs - mean(J_perm)) / sd(J_perm) is reported as a standardized
# effect size (cf. this project's beta_std convention for the regression screens).
#
# CAVEAT -- seasonality is only partly controlled. Within-trap shuffling
# preserves which site/trap/lure a sample belongs to but not WHEN in the
# season it was taken, so a positive hit can still partly reflect both
# partners peaking in spring rather than a direct association. A
# date-conditional version (block on collection date) is the natural
# follow-up if a pair looks interesting.
#
# Read depth: presence/absence uses the RAW ASV table (not rarefied) so the
# rare-tail Ophiostomatales detections are kept. num_reads spans ~100x across
# samples, so ASV-level results in particular should be read with that in
# mind; the genus roll-up is the most robust view.
##############################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(permute)

set.seed(1)

MIN_PREV <- 3       # min samples a fungal species/ASV unit must occur in
N_PERM   <- 1999

out_fig_dir  <- "figures/insect_exploratory/target_genera"
out_data_dir <- "data/insect_exploratory/target_genera"
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("Ips grandicollis", "Ips pini", "Dendroctonus valens")

## ---- 1. Load + match insect and fungal data (same idiom as -----------------
##         compare_insect_fungi/insect_fungal_pairwise_taxon_association.r) ----

sp_tab           <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")
fungal_full      <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
fungal_taxonomy  <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv",
                               row.names = 1, check.names = FALSE)

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
stopifnot(all(targets %in% sp_tab$Finest.ID))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

id_map <- insect_meta_full %>% filter(col_names %in% rownames(sp_tab.t))
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)
cat(length(shared_ids), "samples have both insect and fungal data.\n")

insect <- sp_tab.t[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(insect) == meta$sample_id),
          all(rownames(fungal) == meta$sample_id))

## ---- 2. Beetle + Ophiostomatales presence/absence ------------------------

beetle_pa <- (insect[, targets, drop = FALSE] > 0) * 1L
cat("\nBeetle presence (of", nrow(beetle_pa), "samples):\n")
print(colSums(beetle_pa))

oph_asv <- rownames(fungal_taxonomy)[!is.na(fungal_taxonomy$Order) &
                                       fungal_taxonomy$Order == "o__Ophiostomatales"]
oph_asv <- intersect(oph_asv, colnames(fungal))
oph_pa  <- (fungal[, oph_asv, drop = FALSE] > 0) * 1L
cat("\n", length(oph_asv), "Ophiostomatales ASVs present in the matched samples;",
    "any-Ophiostomatales present in", sum(rowSums(oph_pa) > 0), "samples.\n")

gen <- sub("^g__", "", fungal_taxonomy[oph_asv, "Genus"])
spp <- sub("^s__", "", fungal_taxonomy[oph_asv, "Species"])
gen[is.na(gen)] <- "unclassified"

# Roll a set of ASV columns up to a single presence/absence vector.
rollup <- function(asv_cols) as.integer(rowSums(oph_pa[, asv_cols, drop = FALSE]) > 0)

# (1) species grid -- named species only
sp_label <- ifelse(!is.na(spp) & gen != "unclassified", paste(gen, spp), NA)
species_units <- tapply(seq_along(oph_asv), sp_label, function(ix) ix)
species_pa <- sapply(species_units, rollup)
rownames(species_pa) <- rownames(oph_pa)
species_prev <- colSums(species_pa)
species_pa <- species_pa[, species_prev >= MIN_PREV, drop = FALSE]
cat("\nSpecies grid:", ncol(species_pa), "named Ophiostomatales species at prevalence >=", MIN_PREV, "\n")

# (2) genus grid + "Any Ophiostomatales"
genus_units <- tapply(seq_along(oph_asv), gen, function(ix) ix)
genus_pa <- sapply(genus_units, rollup)
rownames(genus_pa) <- rownames(oph_pa)
genus_pa <- cbind(genus_pa, `Any Ophiostomatales` = as.integer(rowSums(oph_pa) > 0))
colnames(genus_pa) <- sub("^unclassified$", "Ophiostomatales (unclassified genus)", colnames(genus_pa))
cat("Genus grid:", ncol(genus_pa), "units (", paste(colnames(genus_pa), collapse = ", "), ")\n")

# (3) ASV grid
asv_prev <- colSums(oph_pa)
asv_pa <- oph_pa[, asv_prev >= MIN_PREV, drop = FALSE]
asv_meta <- tibble(fungal_unit = colnames(asv_pa),
                   genus  = sub("^g__", "", fungal_taxonomy[colnames(asv_pa), "Genus"]),
                   species = sub("^s__", "", fungal_taxonomy[colnames(asv_pa), "Species"]))
cat("ASV grid:", ncol(asv_pa), "Ophiostomatales ASVs at prevalence >=", MIN_PREV, "\n")

# reference presence/absence matrix
write.csv(as.data.frame(cbind(beetle_pa, genus_pa)) %>% rownames_to_column("sample_id"),
          file.path(out_data_dir, "beetle_ophiostomatales_presence_absence.csv"),
          row.names = FALSE)

## ---- 3. Permutation index set (within-trap shuffle) --------------------------

ctrl <- how(blocks = meta$site,
            plots  = Plots(strata = meta$trap_id, type = "none"),
            within = Within(type = "free"),
            nperm  = N_PERM)
perm_mat <- shuffleSet(nrow(meta), control = ctrl)      # N_PERM x n
stopifnot(nrow(perm_mat) == N_PERM)

## ---- 4. Jaccard + permutation + Fisher for one grid ------------------------
# jaccard_grid_test() lives in cooccurrence_lib.R and is shared with
# scolytinae_ophiostomatales_cooccurrence.R.

source("insect_exploratory/cooccurrence_lib.R")

res_species <- jaccard_grid_test(beetle_pa, species_pa, perm_mat, "species")
res_genus   <- jaccard_grid_test(beetle_pa, genus_pa,   perm_mat, "genus")
res_asv     <- jaccard_grid_test(beetle_pa, asv_pa,     perm_mat, "asv") %>%
  left_join(asv_meta, by = "fungal_unit") %>%
  relocate(genus, species, .after = fungal_unit)

write.csv(res_species, file.path(out_data_dir, "beetle_ophiostomatales_jaccard.species.csv"), row.names = FALSE)
write.csv(res_genus,   file.path(out_data_dir, "beetle_ophiostomatales_jaccard.genus.csv"),   row.names = FALSE)
write.csv(res_asv,     file.path(out_data_dir, "beetle_ophiostomatales_jaccard.asv.csv"),     row.names = FALSE)

cat("\n--- Pairs with q_perm < 0.10 ---\n")
hits <- bind_rows(res_species, res_genus, res_asv %>% select(any_of(names(res_species)))) %>%
  filter(q_perm < 0.10) %>% arrange(q_perm)
if (nrow(hits) == 0) cat("(none)\n") else print(as.data.frame(hits), digits = 3)

cat("\n--- Strongest raw co-occurrences (species grid, by Jaccard) ---\n")
print(as.data.frame(res_species %>% arrange(desc(jaccard)) %>% head(10)), digits = 3)

## ---- 5. Heatmaps --------------------------------------------------------------

# Cell colour is the standardized effect (observed Jaccard vs the within-trap
# permutation null) so the eye is not just tracking prevalence; the printed
# number is the Jaccard value itself, with the shared-presence count beneath.
heat <- function(res, title, ylab_col = "fungal_unit", min_total_shared = 2) {
  d <- res %>%
    mutate(fu = .data[[ylab_col]],
           star  = ifelse(q_perm < 0.10, "*", ""),
           label = ifelse(shared_n > 0,
                          sprintf("%.2f%s (%d)", jaccard, star, shared_n), ""))
  # drop fungal units with almost no shared presence with any beetle -- they
  # carry no co-occurrence information and just add empty rows (full grid is
  # kept in the CSV).
  keep_fu <- d %>% group_by(fu) %>% summarise(s = sum(shared_n)) %>%
    filter(s >= min_total_shared) %>% pull(fu)
  d <- d %>% filter(fu %in% keep_fu)
  # order fungal units by total shared presence across the three beetles
  ord <- d %>% group_by(fu) %>% summarise(s = sum(shared_n)) %>% arrange(s) %>% pull(fu)
  d <- d %>% mutate(fu = factor(fu, levels = ord),
                    beetle = factor(beetle, levels = targets))
  lim <- max(abs(d$z_score), na.rm = TRUE)
  ggplot(d, aes(beetle, fu, fill = z_score)) +
    geom_tile(colour = "grey85") +
    geom_tile(data = subset(d, q_perm < 0.10),
              colour = "black", linewidth = 0.8, fill = NA) +
    geom_text(aes(label = label), size = 2.6) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey88",
                         name = "Jaccard vs null\n(z-score)") +
    scale_x_discrete(labels = function(x) sub(" ", "\n", x)) +
    labs(title = title, subtitle = sub_txt, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 8),
          plot.subtitle = element_text(size = 7.5, colour = "grey30"))
}

sub_txt <- paste(
  "fill = standardized effect (obs Jaccard vs within-trap permutation null);",
  "\ntext = Jaccard (shared-presence samples); black box + * = q_perm < 0.10,",
  N_PERM, "permutations")

p_sp <- heat(res_species, "Target beetle species x Ophiostomatales species")
p_gn <- heat(res_genus,   "Target beetle species x Ophiostomatales genus")
p_av <- heat(res_asv %>% mutate(lab = paste0(fungal_unit, "  ",
              ifelse(is.na(species), ifelse(is.na(genus), "", paste0(genus, " sp.")),
                     paste(genus, species)))),
             "Target beetle species x Ophiostomatales ASV", ylab_col = "lab",
             min_total_shared = 4)

ggsave(file.path(out_fig_dir, "beetle_ophiostomatales_jaccard.species.png"), p_sp,
       width = 7, height = 5, dpi = 150)
ggsave(file.path(out_fig_dir, "beetle_ophiostomatales_jaccard.genus.png"), p_gn,
       width = 7, height = 4.3, dpi = 150)
ggsave(file.path(out_fig_dir, "beetle_ophiostomatales_jaccard.asv.png"), p_av,
       width = 8, height = 5.5, dpi = 150)
dims <- list(species = c(7, 5), genus = c(7, 4.3), asv = c(8, 5.5))
for (nm in names(dims)) {
  pp <- get(c(species = "p_sp", genus = "p_gn", asv = "p_av")[nm])
  ggsave(file.path(out_fig_dir, paste0("beetle_ophiostomatales_jaccard.", nm, ".pdf")),
         pp, width = dims[[nm]][1], height = dims[[nm]][2])
}

cat("\nDone. Outputs in", out_data_dir, "and", out_fig_dir, "\n")
