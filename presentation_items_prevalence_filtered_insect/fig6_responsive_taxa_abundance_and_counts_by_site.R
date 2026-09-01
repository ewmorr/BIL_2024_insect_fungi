##############################################################################
# Presentation figure 6 -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of presentation_items/fig6_responsive_taxa_abundance_
# and_counts_by_site.R. Same 2x2 layout and same underlying question -- how
# much of the fungal community's sequence abundance and richness sits in the
# "responsive" (date- or insect-associated) category, by site -- rebuilt on
# the Lineage-B matched dataset and, per user request, on a COMBINED
# significance definition for each association type instead of Lineage A's
# single-test definition:
#
#   "Date-associated"   = q<0.10 on EITHER the linear OR the quadratic date
#                          term (data/2024_fungi/fungal_taxa_date_
#                          association.all_taxa.csv, lineage-invariant --
#                          same file Lineage A's fig6 reads, but Lineage A
#                          only used the linear term; the quadratic term
#                          didn't exist yet when that figure was built).
#   "Insect-associated"  = q<0.10 on EITHER insect_PCoA1 OR insect_PCoA2,
#                          factor(date)-adjusted (fungal_insect_association_
#                          results.PCoA{1,2}_plus_factordate.csv -- the same
#                          deseasonalized files fig3/fig4/fig5 use, in place
#                          of Lineage A's single permissive no-lure PCoA1
#                          test). Both files test the SAME top-200 CoCA-
#                          loading candidates, and their significant hits
#                          don't overlap (17 PCoA1 + 9 PCoA2 = 26 union), so
#                          this is a straightforward union, not a
#                          double-counted overlap.
#
# "Not associated" is everything else, same as Lineage A: taxa tested-and-
# not-significant on every relevant term, and taxa never tested (below the
# >=5-sample prevalence filter, or outside the top-200 CoCA candidates for
# the insect screens).
#
# Panel a: date-associated vs. not, relative sequence abundance, by site
# Panel b: insect-associated vs. not, relative sequence abundance, by site
# Panel c: date-associated vs. not, ASV count (richness), by site
# Panel d: insect-associated vs. not, ASV count (richness), by site
#
# Matched dataset construction mirrors fig2's Lineage-B port exactly (all
# families, individual-taxon >=5-sample prevalence filter, Kingdom==Fungi
# fungal ASVs, samples shared between both tables) -- see that script for the
# full rationale. Site labels use the raw Site column (no "Pease" ->
# "Pease Airport" rename), matching fig2's Lineage-B convention rather than
# Lineage A's fig6.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)

out_fig_dir <- "figures/presentation_items_prevalence_filtered_insect"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

lb_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"

q_threshold <- 0.10

## ---- Load + match the Lineage-B insect/fungal dataset ----------------------
# Same construction as fig2_nmds_procrustes.R's Lineage-B port: all families,
# individual-taxon >=5-sample prevalence filter (not Curculionidae+Latridiidae).

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
id_map <- id_map %>% filter(SequenceID %in% shared_ids)
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
fungal <- fungal[, is_fungus, drop = FALSE]

cat("Matched dataset:", nrow(fungal), "samples,", ncol(fungal), "Kingdom==Fungi ASVs.\n")

site_levels <- c("Durham", "Pease", "Manchester Cedar Swamp", "Manchester Airport")
meta <- id_map %>%
  transmute(SequenceID, Site = factor(Site, levels = site_levels))

## ---- The two association screens, each a COMBINED (union) definition ------

date_assoc <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
date_sig <- date_assoc %>% filter(q_value < q_threshold | q_value_quad < q_threshold) %>% pull(taxon)

pcoa1_sig <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv")) %>%
  filter(q_value < q_threshold) %>% pull(taxon)
pcoa2_sig <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv")) %>%
  filter(q_value < q_threshold) %>% pull(taxon)
insect_sig <- union(pcoa1_sig, pcoa2_sig)

cat(length(date_sig), "date-associated (linear or quadratic, union),",
    length(insect_sig), "insect-associated (PCoA1 or PCoA2, union) taxa,",
    "out of", ncol(fungal), "Kingdom==Fungi ASVs.\n")

## ---- Taxon-level long table (one row per sample x present taxon) ----------

taxon_long <- as.data.frame(fungal) %>%
  tibble::rownames_to_column("SequenceID") %>%
  pivot_longer(cols = -SequenceID, names_to = "taxon", values_to = "count") %>%
  filter(count > 0) %>%
  left_join(meta, by = "SequenceID")

## ---- Shared plotting logic --------------------------------------------------

pos_color <- "#0072B2"
neg_color <- "grey70"

panel_theme <- function(show_x_text = TRUE, show_legend = TRUE, show_y = TRUE) {
  theme_bw() +
    theme(
      axis.text = element_text(size = 10, color = "black"),
      axis.text.x = if (show_x_text) element_text(angle = 30, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x_text) element_line() else element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = if (show_y) element_text(size = 10, color = "black") else element_blank(),
      axis.title.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      plot.subtitle = element_text(size = 8.5, color = "grey30"),
      plot.tag = element_text(size = 14, face = "bold"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
}

make_abundance_panel <- function(sig_taxa, pos_label, neg_label, title, subtitle, tag, show_y = TRUE) {
  prop <- taxon_long %>%
    mutate(group = ifelse(taxon %in% sig_taxa, pos_label, neg_label)) %>%
    group_by(SequenceID, group) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(SequenceID) %>%
    mutate(sample_total = sum(count), prop = count / sample_total) %>%
    ungroup() %>%
    left_join(meta, by = "SequenceID")

  site_mean <- prop %>%
    group_by(Site, group) %>%
    summarise(mean_prop = mean(prop), .groups = "drop") %>%
    mutate(group = factor(group, levels = c(pos_label, neg_label)))

  ggplot(site_mean, aes(x = Site, y = mean_prop, fill = group)) +
    geom_col(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label))) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(y = "Mean relative sequence abundance", title = title, subtitle = subtitle, tag = tag) +
    panel_theme(show_x_text = FALSE, show_legend = FALSE, show_y = show_y)
}

make_count_panel <- function(sig_taxa, pos_label, neg_label, title, subtitle, tag, show_y = TRUE) {
  site_counts <- taxon_long %>%
    distinct(Site, taxon) %>%
    mutate(group = ifelse(taxon %in% sig_taxa, pos_label, neg_label)) %>%
    count(Site, group) %>%
    mutate(group = factor(group, levels = c(pos_label, neg_label)))

  ggplot(site_counts, aes(x = Site, y = n, fill = group)) +
    geom_col(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y = "ASVs present (count)", title = title, subtitle = subtitle, tag = tag) +
    panel_theme(show_x_text = TRUE, show_legend = TRUE, show_y = show_y)
}

## ---- Build the four panels ---------------------------------------------------

date_subtitle <- paste0(length(date_sig), "/", ncol(fungal), " taxa, q<0.10 linear or quadratic date term")
insect_subtitle <- paste0(length(insect_sig), "/", ncol(fungal), " taxa, q<0.10 PCoA1 or PCoA2, factor(date)-adjusted")

panel_a <- make_abundance_panel(date_sig, "Date-associated", "Not date-associated",
                                 "Date-associated: relative abundance", date_subtitle, "a", show_y = TRUE)
panel_b <- make_abundance_panel(insect_sig, "Insect-associated", "Not insect-associated",
                                 "Insect-associated: relative abundance", insect_subtitle, "b", show_y = FALSE)
panel_c <- make_count_panel(date_sig, "Date-associated", "Not date-associated",
                             "Date-associated: ASV count", date_subtitle, "c", show_y = TRUE)
panel_d <- make_count_panel(insect_sig, "Insect-associated", "Not insect-associated",
                             "Insect-associated: ASV count", insect_subtitle, "d", show_y = FALSE)

## ---- Combine + save -----------------------------------------------------------
# All four panels share a flat (unfaceted) gtable structure, so gtable::
# rbind/cbind aligns axes cleanly -- same approach as Lineage A's fig6.

ga <- ggplotGrob(panel_a); gb <- ggplotGrob(panel_b)
gc <- ggplotGrob(panel_c); gd <- ggplotGrob(panel_d)
left_col <- rbind(ga, gc, size = "max")
right_col <- rbind(gb, gd, size = "max")
combined <- cbind(left_col, right_col, size = "max")

ggsave(file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.png"), combined,
       width = 9, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.pdf"), combined,
       width = 9, height = 8.5, bg = "white")

cat("\nDone. Wrote", file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.png"), "and .pdf\n")
