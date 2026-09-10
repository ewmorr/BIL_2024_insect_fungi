##############################################################################
# Date-conditional co-occurrence: top-scoring Scolytinae x Ophiostomatales
#
# scolytinae_ophiostomatales_cooccurrence.R shuffles beetle presence/absence
# WITHIN TRAP. That controls site/trap/lure but NOT collection date: only the
# beetle vector is permuted while the fungal vector keeps its full seasonal
# signal, so two partners that both simply peak in spring score as
# "co-occurring". The genus heat map from that script showed a whole block of
# early-season Scolytinae turning uniformly positive against nearly every
# Ophiostomatales genus -- the signature of a shared seasonal driver rather
# than specific associations.
#
# This script re-tests that block with the SAME observed Jaccard but a
# date-conditional null: the beetle P/A vector is permuted only AMONG SAMPLES
# THAT SHARE A COLLECTION DATE, so each date's beetle prevalence is held
# fixed and shared seasonality cannot generate signal. Two stratifications:
#
#   within_date         blocks = date            -- season held fixed
#   within_site_x_date  blocks = site x date     -- season AND site held fixed
#                                                  (permutes beetle among the
#                                                  3 lure-traps of a site-date;
#                                                  strict, lower power)
#
# The original within_trap null is recomputed here too so the three sit
# side by side. A pair whose z stays positive under within_date is a genuine
# date-independent candidate; a pair whose z collapses toward 0 was seasonal.
#
# This is the presence/absence analogue of the project's
# insect_fungal_pairwise_taxon_association.residualized.r (which residualizes
# the CLR insect predictor on site + factor(date)); stratified permutation is
# used instead of residualization because Jaccard is a set metric, not a
# regression coefficient.
#
# Shared engine: insect_exploratory/cooccurrence_lib.R (jaccard_grid_test).
##############################################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(permute)

set.seed(1)

MIN_PREV <- 3
N_PERM   <- 1999

out_fig_dir  <- "figures/insect_exploratory/scolytinae"
out_data_dir <- "data/insect_exploratory/scolytinae"
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("Ips grandicollis", "Ips pini", "Dendroctonus valens")

# Top-scoring block from scolytinae_ophiostomatales_cooccurrence.R
# (mean z across the 9 Ophiostomatales genera, within-trap null, >= 1.5):
BLOCK <- c("Pityogenes hopkinsi", "Xyleborinus attenuatus", "Cyclorhipidion pelliculosum",
           "Dendroctonus valens", "Heteroborips seriatus", "Xyloterinus politus",
           "Hylesinus aculeatus", "Hylastes opacus")
# carried through as contrast rows, not part of the block hypothesis set:
REFERENCE <- c("Ips grandicollis", "Ips pini", "Dryocoetes autographus")

## ---- 1. Load + match (same idiom as the non-residualized script) --------

sp_tab           <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")
fungal_full      <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
fungal_taxonomy  <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv",
                               row.names = 1, check.names = FALSE)

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

id_map <- insect_meta_full %>% filter(col_names %in% rownames(sp_tab.t))
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- sp_tab.t[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(insect) == meta$sample_id),
          all(rownames(fungal) == meta$sample_id))
cat(nrow(meta), "matched samples;", n_distinct(meta$date), "collection dates,",
    "date-block sizes:", paste(sort(unique(table(meta$date))), collapse = "/"), "\n")

## ---- 2. Beetle block + Ophiostomatales presence/absence ----------------

beetle_set <- c(BLOCK, REFERENCE)
stopifnot(all(beetle_set %in% colnames(insect)))
beetle_pa <- (insect[, beetle_set, drop = FALSE] > 0) * 1L
cat("\nBeetle presence (of 69):\n"); print(colSums(beetle_pa))

oph_asv <- rownames(fungal_taxonomy)[!is.na(fungal_taxonomy$Order) &
                                       fungal_taxonomy$Order == "o__Ophiostomatales"]
oph_asv <- intersect(oph_asv, colnames(fungal))
oph_pa  <- (fungal[, oph_asv, drop = FALSE] > 0) * 1L

gen <- sub("^g__", "", fungal_taxonomy[oph_asv, "Genus"]); gen[is.na(gen)] <- "unclassified"
spp <- sub("^s__", "", fungal_taxonomy[oph_asv, "Species"])
rollup <- function(ix) as.integer(rowSums(oph_pa[, ix, drop = FALSE]) > 0)

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
                   asv_genus   = sub("^g__", "", fungal_taxonomy[colnames(asv_pa), "Genus"]),
                   asv_species = sub("^s__", "", fungal_taxonomy[colnames(asv_pa), "Species"]))

cat("Fungal grids: genus", ncol(genus_pa), "| species", ncol(species_pa),
    "| ASV", ncol(asv_pa), "\n")

## ---- 3. Three permutation schemes -------------------------------------

schemes <- list(
  within_trap = how(blocks = meta$site,
                    plots  = Plots(strata = meta$trap_id, type = "none"),
                    within = Within(type = "free"), nperm = N_PERM),
  within_date = how(blocks = factor(meta$date),
                    within = Within(type = "free"), nperm = N_PERM),
  within_site_x_date = how(blocks = factor(paste(meta$site, meta$date)),
                           within = Within(type = "free"), nperm = N_PERM)
)
perm_mats <- lapply(schemes, function(ct) shuffleSet(nrow(meta), control = ct))
for (nm in names(perm_mats))
  cat(sprintf("  %-18s %d unique permutations sampled\n", nm,
              nrow(unique(perm_mats[[nm]]))))

## ---- 4. Run the shared engine for every (scheme x grid) --------------

source("insect_exploratory/cooccurrence_lib.R")

grids <- list(genus = genus_pa, species = species_pa, asv = asv_pa)
res <- list()
for (sc in names(perm_mats)) for (g in names(grids)) {
  r <- jaccard_grid_test(beetle_pa, grids[[g]], perm_mats[[sc]], g) %>%
    mutate(scheme = sc, in_block = beetle %in% BLOCK)
  res[[paste(sc, g, sep = ".")]] <- r
}
res_all <- bind_rows(res) %>%
  left_join(asv_meta, by = "fungal_unit") %>%
  mutate(asv_lab = ifelse(grid == "asv",
                          paste0(fungal_unit, "  ",
                                 ifelse(is.na(asv_species),
                                        ifelse(is.na(asv_genus), "", paste0(asv_genus, " sp.")),
                                        paste(asv_genus, asv_species))),
                          fungal_unit))

## ---- 5. Wide comparison table (jaccard fixed; null varies) -----------

cmp <- res_all %>%
  select(grid, beetle, fungal_unit, asv_genus, asv_species, in_block,
         beetle_n, fungus_n, shared_n, jaccard,
         scheme, exp_jaccard_null, z_score, p_pos, p_perm, q_perm) %>%
  pivot_wider(names_from = scheme,
              values_from = c(exp_jaccard_null, z_score, p_pos, p_perm, q_perm)) %>%
  arrange(grid, match(beetle, c(BLOCK, REFERENCE)), desc(jaccard))

write.csv(cmp, file.path(out_data_dir,
          "scolytinae_ophiostomatales_jaccard.date_conditional_comparison.csv"),
          row.names = FALSE)
write.csv(res_all, file.path(out_data_dir,
          "scolytinae_ophiostomatales_jaccard.date_conditional_long.csv"),
          row.names = FALSE)

## ---- 6. Console summary ---------------------------------------------

cat("\n=== Mean z per block beetle (genus grid, 'Any Ophiostomatales' excluded) ===\n")
print(as.data.frame(
  res_all %>% filter(grid == "genus", fungal_unit != "Any Ophiostomatales") %>%
    group_by(beetle, in_block) %>%
    summarise(z_within_trap = mean(z_score[scheme == "within_trap"], na.rm = TRUE),
              z_within_date = mean(z_score[scheme == "within_date"], na.rm = TRUE),
              z_site_x_date = mean(z_score[scheme == "within_site_x_date"], na.rm = TRUE),
              n_pos_p05_date = sum(p_pos[scheme == "within_date"] < 0.05, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(z_within_trap))), digits = 3)

cat("\n=== within_date hits per grid ===\n")
print(as.data.frame(
  res_all %>% filter(scheme == "within_date") %>%
    group_by(grid) %>%
    summarise(n_pairs = n(),
              p05 = sum(p_perm < 0.05, na.rm = TRUE),
              q10 = sum(q_perm < 0.10, na.rm = TRUE), .groups = "drop")))

cat("\n=== Pairs q_perm < 0.10 under within_date (all grids) ===\n")
keep <- res_all %>% filter(scheme == "within_date", p_perm < 0.05) %>%
  arrange(grid, p_perm) %>%
  select(grid, beetle, fungal_unit, in_block, shared_n, jaccard, z_score, p_perm, q_perm, direction)
kq <- keep %>% filter(q_perm < 0.10)
if (nrow(kq) == 0) cat("(none)\n") else print(as.data.frame(kq), digits = 3)

cat("\n=== All pairs p_perm < 0.05 under within_date (all grids) ===\n")
if (nrow(keep) == 0) cat("(none)\n") else print(as.data.frame(keep), digits = 3)

cat("\n=== ... of those, which also reach within_site_x_date p_perm < 0.05 ===\n")
strict <- res_all %>% filter(scheme == "within_site_x_date", p_perm < 0.05) %>%
  semi_join(keep, by = c("grid", "beetle", "fungal_unit")) %>%
  arrange(grid, p_perm) %>%
  select(grid, beetle, fungal_unit, jaccard, z_score, p_perm, q_perm)
if (nrow(strict) == 0) cat("(none)\n") else print(as.data.frame(strict), digits = 3)

## ---- 7. Figures: z under the three nulls, one per fungal grid -------

scheme_lab <- c(within_trap = "within trap\n(site/trap/lure fixed; season NOT)",
                within_date = "within date\n(season fixed)",
                within_site_x_date = "within site x date\n(season + site fixed)")

beetle_lvls <- rev(c(BLOCK, REFERENCE))
beetle_face <- ifelse(beetle_lvls %in% BLOCK, "bold", "italic")

# orient = "fungal_y": fungal unit on y, beetle on x (used for the ASV grid,
# which has too many fungal units to sit on the x axis); otherwise fungal
# unit on x, beetle on y (genus, species).
panel_fig <- function(grid_name, title, fu_col = "fungal_unit", orient = "beetle_y") {
  # beetle axis: block first (best-scoring at the readable end), reference last.
  # on a y axis ggplot draws bottom-up, so reverse there.
  blvls <- if (orient == "fungal_y") c(BLOCK, REFERENCE) else beetle_lvls
  d <- res_all %>% filter(grid == grid_name) %>%
    mutate(fu     = .data[[fu_col]],
           beetle = factor(beetle, levels = blvls),
           scheme = factor(scheme_lab[scheme], levels = scheme_lab),
           star   = ifelse(q_perm < 0.10, "*", ifelse(p_perm < 0.05, "·", "")),
           lab    = ifelse(shared_n >= 3, sprintf("%.2f%s", jaccard, star), star))
  fu_ord <- d %>% group_by(fu) %>% summarise(s = sum(shared_n)) %>%
    arrange(desc(s)) %>% pull(fu)
  d <- d %>% mutate(fu = factor(fu, levels = if (orient == "fungal_y") rev(fu_ord) else fu_ord))
  lim <- max(abs(d$z_score), na.rm = TRUE)
  mapping <- if (orient == "fungal_y") aes(beetle, fu, fill = z_score)
             else                      aes(fu, beetle, fill = z_score)
  gg <- ggplot(d, mapping) +
    geom_tile(colour = "grey88") +
    geom_tile(data = subset(d, q_perm < 0.10), colour = "black", linewidth = 0.6, fill = NA) +
    geom_text(aes(label = lab), size = 2.1) +
    facet_wrap(~ scheme, nrow = 1) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                         midpoint = 0, limits = c(-lim, lim), na.value = "grey90",
                         name = "Jaccard vs\nnull (z)") +
    labs(x = NULL, y = NULL, title = title,
         subtitle = paste("Same observed Jaccard in every panel; only the null changes.",
                          "Number = Jaccard (>=3 shared samples); * q_perm<0.10, dot p_perm<0.05.",
                          "Bold beetle = block; italic = reference.")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.subtitle = element_text(size = 7.5, colour = "grey30"),
          strip.text = element_text(size = 8))
  if (orient == "fungal_y")
    gg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                          face = ifelse(blvls %in% BLOCK, "bold", "italic")),
               axis.text.y = element_text(size = 7))
  else
    gg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
               axis.text.y = element_text(size = 8, face = beetle_face))
}

p_gn <- panel_fig("genus",   "Date-conditional co-occurrence: Scolytinae block x Ophiostomatales genus")
p_sp <- panel_fig("species", "Date-conditional co-occurrence: Scolytinae block x Ophiostomatales species")
p_av <- panel_fig("asv",     "Date-conditional co-occurrence: Scolytinae block x Ophiostomatales ASV",
                  fu_col = "asv_lab", orient = "fungal_y")

ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.genus.png"),   p_gn, width = 15,  height = 5.5, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.species.png"), p_sp, width = 17,  height = 5.5, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.asv.png"),     p_av, width = 12,  height = 9,   dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.genus.pdf"),   p_gn, width = 15,  height = 5.5)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.species.pdf"), p_sp, width = 17,  height = 5.5)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.asv.pdf"),     p_av, width = 12,  height = 9)

## ---- 8. Figure: z shrinkage across the three nulls -------------------

slope <- res_all %>%
  filter(grid == "genus", fungal_unit != "Any Ophiostomatales") %>%
  mutate(scheme = factor(scheme, levels = names(scheme_lab), labels = c("trap", "date", "site x date")))

p2 <- ggplot(slope, aes(scheme, z_score, group = interaction(beetle, fungal_unit),
                        colour = beetle %in% BLOCK)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(alpha = 0.35) + geom_point(size = 0.8, alpha = 0.5) +
  scale_colour_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "block beetle", `FALSE` = "reference"),
                      name = NULL) +
  labs(x = "permutation null (season control increasing ->)", y = "Jaccard vs null (z)",
       title = "Each beetle x genus pair: standardized co-occurrence under the three nulls") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())

ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.shrinkage.png"),
       p2, width = 6.5, height = 5, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_ophiostomatales_jaccard.date_conditional.shrinkage.pdf"),
       p2, width = 6.5, height = 5)

cat("\nDone. Outputs in", out_data_dir, "and", out_fig_dir, "\n")
