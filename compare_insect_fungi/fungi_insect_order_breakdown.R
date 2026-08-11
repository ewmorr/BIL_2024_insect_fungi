# Exploratory: Order-level taxonomic breakdown of the fungal taxa
# significantly (q<0.10) associated with insect_PCoA1, no lure/date
# covariates -- data/compare_insects_fungi_top3axes/fungal_insect_
# association_results.PCoA1_only.no_lure.csv (same 77 hits as fig3/fig4
# panel c, top-200 CoCA-loading candidates). Family was already broken down
# in fig4 (presentation_items/fig4_association_taxonomic_breakdown.R);
# Order is one level coarser -- fewer, broader groups, useful for a quick
# look at whether the pattern holds up (or looks different) at that level.

library(dplyr)
library(ggplot2)

fig_dir <- "figures"

q_threshold <- 0.10

tax <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1) %>%
  tibble::rownames_to_column("taxon") %>%
  transmute(taxon, Order = sub("^o__", "", Order))

res <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(tax, by = "taxon") %>%
  mutate(
    Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order),
    direction = ifelse(t_stat > 0, "PCoA1+ (later-season insect community)", "PCoA1- (bark-beetle-dominated insect community)")
  )

cat(nrow(res), "significant taxa across", n_distinct(res$Order), "orders\n")

counts <- res %>%
  count(Order, direction) %>%
  mutate(signed_n = ifelse(grepl("\\+", direction), n, -n))

order_rank <- counts %>% group_by(Order) %>% summarize(total = sum(n)) %>% arrange(total)
counts <- counts %>% mutate(Order = factor(Order, levels = order_rank$Order))

p <- ggplot(counts, aes(x = Order, y = signed_n, fill = direction)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 3) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = c("PCoA1+ (later-season insect community)" = "#0072B2",
                                "PCoA1- (bark-beetle-dominated insect community)" = "#D55E00")) +
  labs(x = NULL, y = "Number of significant taxa (q<0.10)", fill = "Direction",
       title = "Insect-associated fungal taxa (PCoA1, no lure/date) by Order",
       subtitle = paste0(nrow(res), " significant taxa across ", n_distinct(res$Order), " orders")) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

print(p)
ggsave(file.path(fig_dir, "fungi_insect_order_breakdown.png"), p, width = 8, height = 7.5, dpi = 300, bg = "white")

cat("Wrote", file.path(fig_dir, "fungi_insect_order_breakdown.png"), "\n")
