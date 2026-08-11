# Exploratory: Order-level taxonomic breakdown of the fungal taxa
# significantly (q<0.10) associated with collection date -- data/2024_fungi/
# fungal_taxa_date_association.all_taxa.csv (unbiased screen, every
# prevalence-filtered taxon, no pre-selection; same 2116 hits as fig3/fig4
# panel b). Family was already broken down in fig4 (presentation_items/
# fig4_association_taxonomic_breakdown.R); Order is one level coarser.
# Companion to fungi_insect_order_breakdown.R (same breakdown for the
# insect_PCoA1-associated taxa).
#
# 2116 hits span 102 orders -- too many to show unfolded (unlike the PCoA1
# screen's 31), so only the top n_top orders by hit count are shown
# individually and the rest folded into "Other".

library(dplyr)
library(ggplot2)

fig_dir <- "figures"

q_threshold <- 0.10
n_top <- 15   # + "Other" = 16 categories

tax <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1) %>%
  tibble::rownames_to_column("taxon") %>%
  transmute(taxon, Order = sub("^o__", "", Order))

res <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(tax, by = "taxon") %>%
  mutate(
    Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order),
    direction = ifelse(t_stat > 0, "Later season (t>0)", "Earlier season (t<0)")
  )

cat(nrow(res), "significant taxa across", n_distinct(res$Order), "orders\n")

order_rank <- res %>% count(Order, name = "total") %>% arrange(desc(total))
top_orders <- head(order_rank$Order, n_top)
res <- res %>% mutate(Order = ifelse(Order %in% top_orders, Order, "Other"))
order_levels <- c(top_orders, "Other")

counts <- res %>%
  count(Order, direction) %>%
  mutate(
    Order = factor(Order, levels = order_levels),
    signed_n = ifelse(direction == "Later season (t>0)", n, -n)
  )

p <- ggplot(counts, aes(x = Order, y = signed_n, fill = direction)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 3) +
  scale_x_discrete(limits = rev(order_levels)) +
  coord_flip(clip = "off") +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.08, 0.08))) +
  scale_fill_manual(values = c("Later season (t>0)" = "#0072B2", "Earlier season (t<0)" = "#D55E00")) +
  labs(x = NULL, y = "Number of significant taxa (q<0.10)", fill = "Direction",
       title = "Date-associated fungal taxa by Order",
       subtitle = paste0(nrow(res), " significant taxa (unbiased screen) across ",
                          n_distinct(order_rank$Order), " orders; top ", n_top, " shown")) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

print(p)
ggsave(file.path(fig_dir, "fungi_date_order_breakdown.png"), p, width = 8, height = 7.5, dpi = 300, bg = "white")

cat("Wrote", file.path(fig_dir, "fungi_date_order_breakdown.png"), "\n")
