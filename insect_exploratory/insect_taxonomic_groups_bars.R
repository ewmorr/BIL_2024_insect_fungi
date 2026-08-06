library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

sp_tab = read.csv("data/2024_insect_data/insect_species_tab.csv")
head(sp_tab)
str(sp_tab)

# some data ops
# first get sum on orders
sp_tab %>% 
    dplyr::group_by(Order) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -Order, names_to = "sampleID", values_to = "Count") %>%
    group_by(Order) %>%
    summarize(Count = sum(Count)) -> global_order_sum

    global_order_sum %>%
    ggplot(., aes(x = reorder(Order, Count), y = Count)) +
        geom_col() +
        theme_bw() +
        scale_y_log10() +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1),
            axis.title.x = element_blank()
        ) -> p1
p1
pdf("figures/insect_comps_2024/order_global.pdf", height = 3.5)
p1
dev.off()

sp_tab %>% 
    filter(Order == "Coleoptera") %>%
    dplyr::group_by(Family) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -Family, names_to = "sampleID", values_to = "Count") %>%
    group_by(Family) %>%
    summarize(Count = sum(Count)) %>%
    ggplot(., aes(x = reorder(Family, Count), y = Count)) +
        geom_col() +
        theme_bw() +
        scale_y_log10() +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1),
            axis.title.x = element_blank()
        ) -> p2
pdf("figures/insect_comps_2024/Coleoptera_family_global.pdf", height = 3.5, width = 12)
p2
dev.off()

sp_tab %>% 
    filter(Order == "Coleoptera") %>%
    dplyr::group_by(Family) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -Family, names_to = "sampleID", values_to = "Count") %>%
    group_by(Family) %>%
    summarize(Count = sum(Count)) %>%
    filter(Count > 100) %>%
    ggplot(., aes(x = reorder(Family, Count), y = Count)) +
        geom_col() +
        theme_bw() +
        #scale_y_log10() +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1),
            axis.title.x = element_blank()
        ) -> p3

pdf("figures/insect_comps_2024/Coleoptera_family_global_gt100count.pdf", height = 3.5, width = 7)
p3
dev.off()


sp_tab %>% 
    filter(Family == "Curculionidae") %>%
    #filter(Family == "Latridiidae") %>%
    dplyr::group_by(Subfamily) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -Subfamily, names_to = "sampleID", values_to = "Count") %>%
    group_by(Subfamily) %>%
    summarize(Count = sum(Count)) %>%
    filter(Count > 100) %>%
    ggplot(., aes(x = reorder(Subfamily, Count), y = Count)) +
        geom_col() +
        theme_bw() +
        #scale_y_log10() +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1),
            axis.title.x = element_blank()
        ) -> p4

pdf("figures/insect_comps_2024/Curculionidae_subfamily_global.pdf", height = 3.5, width = 7)
p4
dev.off()


sp_tab %>% 
    filter(Subfamily == "Scolytinae") %>%
    dplyr::group_by(Genus) %>%
    summarise(across(where(is.numeric), sum)) %>%
    pivot_longer(cols = -Genus, names_to = "sampleID", values_to = "Count") %>%
    group_by(Genus) %>%
    summarize(Count = sum(Count)) %>%
    filter(Count > 100) %>%
    ggplot(., aes(x = reorder(Genus, Count), y = Count)) +
        geom_col() +
        theme_bw() +
        #scale_y_log10() +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1),
            axis.title.x = element_blank()
        ) -> p5

pdf("figures/insect_comps_2024/Scolytinae_genera_global.pdf", height = 3.5, width = 7)
p5
dev.off()
