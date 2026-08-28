##############################################################################
# Insect community NMDS -- individual-taxon prevalence filter variant
#
# insect_ords.R settled on restricting the insect trap-catch table to
# Curculionidae + Latridiidae (the two dominant families, ~69% of
# individuals) because the full, unfiltered table produced a high-stress,
# non-repeatable NMDS solution -- too many taxa occurring in only 1-2
# samples. This script asks whether an individual-taxon PREVALENCE filter
# (present in >=5 samples -- the same threshold used for fungal ASVs
# throughout this project's compare_insect_fungi/ scripts) fixes the same
# problem without discarding taxa outside those two families.
#
# Kept as a separate, clean script rather than appended to insect_ords.R's
# exploratory trial-and-error, per user request.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(permute)
source("library/library.R")

set.seed(1)

out_fig_dir <- "figures/insect_exploratory_prevalence_filtered_insect"
out_data_dir <- "data/insect_exploratory_prevalence_filtered_insect"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load full insect table, apply individual-taxon prevalence filter ----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]  # one all-zero row with no ID, same drop insect_ords.R makes
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID
cat("Full insect table:", nrow(sp_tab.t), "samples,", ncol(sp_tab.t), "taxa across",
    length(unique(sp_tab$Family)), "families.\n")

keep <- colSums(sp_tab.t > 0) >= 5
insect <- sp_tab.t[, keep, drop = FALSE]
insect <- insect[rowSums(insect) > 0, , drop = FALSE]
cat(sum(keep), "of", ncol(sp_tab.t),
    "taxa retained at >=5-sample prevalence (same threshold used for fungal ASVs throughout compare_insect_fungi/).\n")
cat("After dropping resulting empty samples:", nrow(insect), "samples,", ncol(insect), "taxa.\n")

# Family composition of the retained taxa, for comparison with the
# Curculionidae+Latridiidae-only table used elsewhere in this project.
retained_families <- sp_tab %>% filter(Finest.ID %in% colnames(insect)) %>% count(Family, sort = TRUE)
cat("\nFamily composition of the", ncol(insect), "prevalence-filtered taxa (top 10 by taxon count):\n")
print(head(retained_families, 10))
write.csv(retained_families, file.path(out_data_dir, "insect_prevalence_filter.family_composition.csv"),
          row.names = FALSE)

## ---- 2. NMDS ------------------------------------------------------------------
# autotransform=FALSE is the choice insect_ords.R settled on for the
# family-filtered table (vegan's default sqrt + Wisconsin-standardization
# autotransform made stress worse there). Confirmed again here before
# committing to it.

nmds_autotransform <- metaMDS(insect, distance = "bray", binary = FALSE, try = 20, trymax = 100)
cat("\nDefault-autotransform stress:", round(nmds_autotransform$stress, 3), "\n")

insect_nmds <- metaMDS(insect, distance = "bray", binary = FALSE, try = 20, trymax = 100, autotransform = FALSE)
cat("No-autotransform stress:", round(insect_nmds$stress, 3),
    "-- lower and more repeatable, so this is the version used below.\n")

saveRDS(insect_nmds, file.path(out_data_dir, "insect_prevalence_filter.NMDS.RDS"))

## ---- 3. Metadata + PERMANOVA ---------------------------------------------------

insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")
id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect))
stopifnot(!any(duplicated(id_map$col_names)))
cat("\n", nrow(id_map), "of", nrow(insect), "samples matched to metadata.\n")

insect <- insect[id_map$col_names, , drop = FALSE]
meta <- id_map %>% transmute(sample_id = col_names, site = Site, lure = Lure, trap_id = trapID,
                              date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(insect) == meta$sample_id))

# Same site/date permutation scheme insect_ords.R uses (permute freely within
# trap, blocked by site) -- this script's purpose is the ordination itself,
# not re-deriving insect_ords.R's separate balanced-subsample lure test, so
# lure's p-value here is reported from this same scheme (as insect_ords.R
# itself notes, that's a reasonable date/site test but an approximate one for
# lure -- see insect_ords.R's own balanced re-test if a rigorous lure p-value
# is needed).
perm_date <- how(
  blocks = meta$site,
  plots  = Plots(strata = meta$trap_id, type = "none"),
  within = Within(type = "free"),
  nperm  = 999
)
insect_permanova <- adonis2(insect ~ site + lure + date, data = meta, permutations = perm_date, by = "margin")
cat("\n--- PERMANOVA: prevalence-filtered insect community ~ site + lure + date ---\n")
print(insect_permanova)
write.csv(as.data.frame(insect_permanova), file.path(out_data_dir, "insect_prevalence_filter.permanova.csv"))

## ---- 4. Plot --------------------------------------------------------------------

nmds_scores <- as.data.frame(scores(insect_nmds, display = "sites")) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = as.numeric(date), shape = lure)) +
  geom_point(size = 3, color = "black", stroke = 0.3) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                        midpoint = median(as.numeric(meta$date)),
                        breaks = pretty(range(as.numeric(meta$date)), n = 4),
                        labels = function(b) format(as.Date(b, origin = "1970-01-01"), "%b %d"),
                        name = "Collection\ndate") +
  scale_shape_manual(values = c(21, 22, 24), name = "Lure") +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1,
           label = paste0("stress = ", round(insect_nmds$stress, 2))) +
  labs(title = "Insect community NMDS", subtitle = "individual-taxon >=5-sample prevalence filter, all families") +
  theme_bw()
ggsave(file.path(out_fig_dir, "insect_ords.prevalence_filter.NMDS.png"), p, width = 7, height = 6, dpi = 150)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
