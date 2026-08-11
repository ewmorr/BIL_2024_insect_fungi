library(dplyr)

fungi_tax = read.table("data/2024_fungi/ASVs_taxonomy.tsv")

fungi_date = read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
fungi_insect_no_date_no_lure = read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv")
fungi_insect_date_lure = read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_plus_date.csv")
fungi_insect_date_no_lure = read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_plus_date.no_lure.csv")
fungi_insect_lure_no_date = read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.csv")

fungi_tax$taxon = row.names(fungi_tax)


fungi_date %>% 
    filter(q_value < 0.1) %>% 
    arrange(t_stat) %>%
    left_join(., fungi_tax) -> fungi_date
fungi_insect_no_date_no_lure %>% 
    filter(q_value < 0.1) %>% 
    arrange(t_stat) %>%
    left_join(., fungi_tax) -> fungi_insect_no_date_no_lure
fungi_insect_date_lure %>% 
    filter(q_value < 0.1) %>% 
    arrange(t_stat) %>%
    left_join(., fungi_tax) -> fungi_insect_date_lure
fungi_insect_date_no_lure %>% 
    filter(q_value < 0.1) %>% 
    arrange(t_stat) %>%
    left_join(., fungi_tax) -> fungi_insect_date_no_lure
fungi_insect_lure_no_date %>% 
    filter(q_value < 0.1) %>% 
    arrange(t_stat) %>%
    left_join(., fungi_tax) -> fungi_insect_lure_no_date

nrow(fungi_date)
# 2116
nrow(fungi_insect_no_date_no_lure)
# 77
nrow(fungi_insect_date_lure)
# 0
nrow(fungi_insect_date_no_lure)
# 0
nrow(fungi_insect_lure_no_date)
# 63
# partialling out lure removes 14 hits. Are they among the top hits?

sum(fungi_insect_no_date_no_lure$taxon %in% fungi_insect_lure_no_date$taxon)
# 59 of the non-partialled hits are shared in the lure partialled
sum(!fungi_insect_lure_no_date$taxon %in% fungi_insect_no_date_no_lure$taxon)
# 4
fungi_insect_lure_no_date[!fungi_insect_lure_no_date$taxon %in% fungi_insect_no_date_no_lure$taxon,]
plot(fungi_insect_lure_no_date$t_stat)
# the four that are picked up in lure partial but not unpartialled are not particularly strong t-stats (-2.7 - -3.4)
# Candida ponderosae
# 
fungi_insect_no_date_no_lure %>% filter(t_stat > 0)

