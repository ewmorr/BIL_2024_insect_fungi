##############################################################################
# Presence/absence co-occurrence:
#   ALL Scolytinae species  x  Ophiostomatales fungi
#
# The subfamily-wide version of target_species_ophiostomatales_cooccurrence.R
# (which is kept as the focused Ips grandicollis / Ips pini / Dendroctonus
# valens view). Same method exactly:
#
#   - Jaccard similarity J = a / (a + b + c); shared absence excluded.
#   - beetle x fungal-taxon pairs ONLY (no beetle-beetle, no fungus-fungus).
#   - Ophiostomatales rolled up three ways: named species (prev >= MIN_PREV),
#     the 8 genera + "Any Ophiostomatales", and individual ASVs (prev >= MIN_PREV).
#   - restricted permutation null: beetle P/A shuffled N_PERM times WITHIN trap
#     (permute::how, blocks = site, plots = trap, within = free); Jaccard
#     recomputed each time; p_perm = 2*min(p_pos, p_neg), BH within grid.
#   - Fisher's exact on the 2x2 as a structure-blind cross-check.
#   - shared engine: insect_exploratory/cooccurrence_lib.R
#
# Beetle set: every Scolytinae species present in >= BEETLE_MIN_PREV of the 69
# insect+fungal samples (the three original targets are always included).
#
# CAVEAT (unchanged): within-trap shuffling controls site/trap/lure but not
# WHEN in the season a sample was taken. Both Scolytinae and Ophiostomatales
# are spring-weighted, so a positive hit can partly reflect shared phenology
# rather than a direct association. Date-conditional follow-up = residualize /
# block on collection date (cf. insect_fungal_pairwise_taxon_association.residualized.r).
##############################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(permute)

set.seed(1)

MIN_PREV        <- 3   # min samples a fungal species/ASV unit must occur in
BEETLE_MIN_PREV <- 3   # min samples a Scolytinae species must occur in
N_PERM          <- 1999

out_fig_dir  <- "figures/insect_exploratory/scolytinae"
out_data_dir <- "data/insect_exploratory/scolytinae"
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("Ips grandicollis", "Ips pini", "Dendroctonus valens")

## ---- 1. Load + match insect and fungal data ------------------------------

sp_tab           <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")
fungal_full      <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
fungal_taxonomy  <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv",
                               row.names = 1, check.names = FALSE)

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
scol_species <- sp_tab %>% filter(Subfamily == "Scolytinae") %>% pull(Finest.ID)
stopifnot(all(targets %in% scol_species))

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

## ---- 2. Beetle (Scolytinae) + Ophiostomatales presence/absence ----------

scol_mat <- insect[, scol_species, drop = FALSE]
scol_prev <- colSums(scol_mat > 0)
keep_beetle <- union(names(scol_prev)[scol_prev >= BEETLE_MIN_PREV], targets)
# order beetle columns by prevalence (desc) for readability
keep_beetle <- names(sort(scol_prev[keep_beetle], decreasing = TRUE))
beetle_pa <- (scol_mat[, keep_beetle, drop = FALSE] > 0) * 1L

cat("\n", ncol(beetle_pa), "of", length(scol_species),
    "Scolytinae species retained at >=", BEETLE_MIN_PREV, "-sample prevalence",
    "(plus any of the 3 targets below that).\n", sep = "")
cat("Beetle prevalence range:", paste(range(colSums(beetle_pa)), collapse = "-"), "\n")

oph_asv <- rownames(fungal_taxonomy)[!is.na(fungal_taxonomy$Order) &
                                       fungal_taxonomy$Order == "o__Ophiostomatales"]
oph_asv <- intersect(oph_asv, colnames(fungal))
oph_pa  <- (fungal[, oph_asv, drop = FALSE] > 0) * 1L
cat(length(oph_asv), "Ophiostomatales ASVs present; any-Ophiostomatales in",
    sum(rowSums(oph_pa) > 0), "samples.\n")

gen <- sub("^g__", "", fungal_taxonomy[oph_asv, "Genus"])
spp <- sub("^s__", "", fungal_taxonomy[oph_asv, "Species"])
gen[is.na(gen)] <- "unclassified"
rollup <- function(asv_cols) as.integer(rowSums(oph_pa[, asv_cols, drop = FALSE]) > 0)

sp_label <- ifelse(!is.na(spp) & gen != "unclassified", paste(gen, spp), NA)
species_pa <- sapply(tapply(seq_along(oph_asv), sp_label, function(ix) ix), rollup)
rownames(species_pa) <- rownames(oph_pa)
species_pa <- species_pa[, colSums(species_pa) >= MIN_PREV, drop = FALSE]

genus_pa <- sapply(tapply(seq_along(oph_asv), gen, function(ix) ix), rollup)
rownames(genus_pa) <- rownames(oph_pa)
genus_pa <- cbind(genus_pa, `Any Ophiostomatales` = as.integer(rowSums(oph_pa) > 0))
colnames(genus_pa) <- sub("^unclassified$", "Ophiostomatales (unclassified genus)", colnames(genus_pa))

asv_pa <- oph_pa[, colSums(oph_pa) >= MIN_PREV, drop = FALSE]
asv_meta <- tibble(fungal_unit = colnames(asv_pa),
                   genus   = sub("^g__", "", fungal_taxonomy[colnames(asv_pa), "Genus"]),
                   species = sub("^s__", "", fungal_taxonomy[colnames(asv_pa), "Species"]))

cat("Fungal grids: species", ncol(species_pa), "| genus", ncol(genus_pa),
    "| ASV", ncol(asv_pa), "\n")

## ---- 3. Permutation index set (within-trap shuffle) --------------------

ctrl <- how(blocks = meta$site,
            plots  = Plots(strata = meta$trap_id, type = "none"),
            within = Within(type = "free"),
            nperm  = N_PERM)
perm_mat <- shuffleSet(nrow(meta), control = ctrl)
stopifnot(nrow(perm_mat) == N_PERM)

## ---- 4. Run the shared Jaccard/permutation/Fisher engine --------------

source("insect_exploratory/cooccurrence_lib.R")

res_species <- jaccard_grid_test(beetle_pa, species_pa, perm_mat, "species")
res_genus   <- jaccard_grid_test(beetle_pa, genus_pa,   perm_mat, "genus")
res_asv     <- jaccard_grid_test(beetle_pa, asv_pa,     perm_mat, "asv") %>%
  left_join(asv_meta, by = "fungal_unit") %>%
  relocate(genus, species, .after = fungal_unit)

add_target_flag <- function(d) d %>% mutate(beetle_is_target = beetle %in% targets)
res_species <- add_target_flag(res_species)
res_genus   <- add_target_flag(res_genus)
res_asv     <- add_target_flag(res_asv)

write.csv(res_species, file.path(out_data_dir, "scolytinae_ophiostomatales_jaccard.species.csv"), row.names = FALSE)
write.csv(res_genus,   file.path(out_data_dir, "scolytinae_ophiostomatales_jaccard.genus.csv"),   row.names = FALSE)
write.csv(res_asv,     file.path(out_data_dir, "scolytinae_ophiostomatales_jaccard.asv.csv"),     row.names = FALSE)
write.csv(as.data.frame(cbind(beetle_pa, genus_pa)) %>% rownames_to_column("sample_id"),
          file.path(out_data_dir, "scolytinae_ophiostomatales_presence_absence.csv"), row.names = FALSE)

all_res <- bind_rows(res_species, res_genus,
                     res_asv %>% select(any_of(names(res_species))))

cat("\n--- Pairs with q_perm < 0.10 (any grid) ---\n")
h <- all_res %>% filter(q_perm < 0.10) %>% arrange(q_perm) %>%
  select(grid, beetle, fungal_unit, beetle_n, fungus_n, shared_n, jaccard,
         z_score, p_perm, q_perm, p_fisher, q_fisher, direction)
if (nrow(h) == 0) cat("(none)\n") else print(as.data.frame(h), digits = 3)

cat("\n--- Pairs with p_perm < 0.05 (uncorrected), genus + species grids ---\n")
print(as.data.frame(
  bind_rows(res_species, res_genus) %>% filter(p_perm < 0.05) %>%
    arrange(p_perm) %>%
    select(grid, beetle, fungal_unit, shared_n, jaccard, z_score, p_perm, q_perm, direction)),
  digits = 3)

cat("\n--- Scolytinae species ranked by mean z-score vs Ophiostomatales genera ---\n")
print(as.data.frame(
  res_genus %>% filter(fungal_unit != "Any Ophiostomatales") %>%
    group_by(beetle, beetle_n, beetle_is_target) %>%
    summarise(mean_z = mean(z_score, na.rm = TRUE),
              n_pos_p05 = sum(p_pos < 0.05, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_z))), digits = 3)

## ---- 5. Heat maps: beetle (rows) x fungal unit (cols) -----------------
# Bigger grid than the 3-species version, so beetles go on the y axis and
# fungal taxa on the x axis. Cell colour = standardized effect (z of observed
# Jaccard vs the within-trap null); printed number = Jaccard for cells with
# >= 3 shared-presence samples; black outline + * = q_perm < 0.10.

sub_txt <- "Beetles ordered by mean z across genera; bold = one of the 3 original targets"
cap_txt <- paste0(
  "Fill = observed Jaccard vs within-trap permutation null (z-score); grey = null s.d. 0.\n",
  "Number = Jaccard for cells with >= 3 shared-presence samples; ",
  "* / black box = q_perm < 0.10 (", N_PERM, " perms, BH within grid).")

heat <- function(res, title, xlab_col = "fungal_unit", min_col_shared = 1) {
  d <- res %>% mutate(fu = .data[[xlab_col]])
  # keep fungal columns with at least some shared presence (CSV keeps all)
  keep_fu <- d %>% group_by(fu) %>% summarise(s = sum(shared_n)) %>%
    filter(s >= min_col_shared) %>% pull(fu)
  d <- d %>% filter(fu %in% keep_fu)
  # beetles ordered by mean z-score (signal), not raw prevalence
  b_ord  <- d %>% group_by(beetle) %>% summarise(mz = mean(z_score, na.rm = TRUE)) %>%
    arrange(mz) %>% pull(beetle)
  fu_ord <- d %>% group_by(fu) %>% summarise(s = sum(shared_n)) %>% arrange(desc(s)) %>% pull(fu)
  tgt <- unique(d$beetle[d$beetle_is_target])
  d <- d %>% mutate(beetle = factor(beetle, levels = b_ord),
                    fu     = factor(fu, levels = fu_ord),
                    star   = ifelse(q_perm < 0.10, "*", ""),
                    lab    = ifelse(shared_n >= 3, sprintf("%.2f%s", jaccard, star), star))
  lim <- max(abs(d$z_score), na.rm = TRUE)
  ggplot(d, aes(fu, beetle, fill = z_score)) +
    geom_tile(colour = "grey88") +
    geom_tile(data = subset(d, q_perm < 0.10), colour = "black", linewidth = 0.7, fill = NA) +
    geom_text(aes(label = lab), size = 2.3) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey90",
                         name = "Jaccard vs\nnull (z)") +
    labs(title = title, subtitle = sub_txt, caption = cap_txt, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8, face = ifelse(b_ord %in% tgt, "bold", "plain")),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          plot.caption = element_text(size = 7, colour = "grey30", hjust = 0),
          plot.caption.position = "plot",
          plot.margin = margin(6, 6, 14, 6),
          plot.title = element_text(size = 12))
}

p_gn <- heat(res_genus,   "All Scolytinae x Ophiostomatales genus")
p_sp <- heat(res_species, "All Scolytinae x Ophiostomatales species")
p_av <- heat(res_asv %>% mutate(lab_fu = paste0(fungal_unit, "  ",
               ifelse(is.na(species), ifelse(is.na(genus), "", paste0(genus, " sp.")),
                      paste(genus, species)))),
             "All Scolytinae x Ophiostomatales ASV", xlab_col = "lab_fu",
             min_col_shared = 5)

ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.genus.png"),   p_gn, width = 8.5, height = 9, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.species.png"), p_sp, width = 11,  height = 9, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.asv.png"),     p_av, width = 12,  height = 9, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.genus.pdf"),   p_gn, width = 8.5, height = 9)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.species.pdf"), p_sp, width = 11,  height = 9)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.asv.pdf"),     p_av, width = 12,  height = 9)

cat("\nDone. Outputs in", out_data_dir, "and", out_fig_dir, "\n")
