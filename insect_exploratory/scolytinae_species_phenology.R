##############################################################################
# Scolytinae species abundance across the season
#
# Extends target_genus_phenology.R (Ips + Dendroctonus only) to EVERY species
# in the subfamily Scolytinae -- the bark/ambrosia beetles, the group the EDRR
# traps are really sampling. The three original target species (Ips
# grandicollis, Ips pini, Dendroctonus valens) are a subset of what is plotted
# here; target_genus_phenology.R is kept as the focused two-genus view.
#
# Outputs:
#   - species x date abundance heat map, faceted by lure (all 45 Scolytinae
#     species; rows ordered by abundance-weighted mean collection date, so the
#     plot reads as an early-to-late-season progression)
#   - subfamily-level totals: individuals and species richness per date x lure
#
# Counts are pooled across the four sites (one trap per site x lure). Data
# source: data/2024_insect_data/EDRRPatho_NHTraps_InsectData_MFD.csv (raw
# long-format trap catch).
#
# Lure codes: E = Ethanol, EA = Alpha-pinene + Ethanol, I = Ips lure.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

out_fig_dir  <- "figures/insect_exploratory/scolytinae"
out_data_dir <- "data/insect_exploratory/scolytinae"
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

targets <- c("Ips grandicollis", "Ips pini", "Dendroctonus valens")

lure_levels <- c("E", "EA", "I")
lure_labels <- c(E  = "Ethanol",
                 EA = "Alpha-pinene + Ethanol",
                 I  = "Ips lure")

## ---- 1. Load + tidy the raw long-format trap catch --------------------------

dat <- read.csv("data/2024_insect_data/EDRRPatho_NHTraps_InsectData_MFD.csv",
                check.names = FALSE)
names(dat) <- make.names(names(dat))
dat <- dat %>% filter(!is.na(Trap.Code), Trap.Code != "")
dat[dat == ""] <- NA
dat <- dat %>%
  mutate(date  = mdy(Collection.Date),
         Count = suppressWarnings(as.numeric(Count)))
n_blank <- sum(is.na(dat$Count))
if (n_blank > 0) {
  cat("Dropping", n_blank, "catch record(s) with no Count entered.\n")
  dat <- dat %>% filter(!is.na(Count))
}
stopifnot(!any(is.na(dat$date)))

## ---- 2. Restrict to Scolytinae -------------------------------------------

scol <- dat %>% filter(Subfamily == "Scolytinae")
stopifnot(all(!is.na(scol$Species)))          # every Scolytinae record is ID'd to species
cat(n_distinct(scol$Species), "Scolytinae species,",
    sum(scol$Count), "individuals, across",
    n_distinct(scol$Genus), "genera.\n")

all_dates <- sort(unique(dat$date))

## ---- 3. Pool across sites: per species x lure x date -----------------------

summ <- scol %>%
  group_by(Genus, Species, Lure, date) %>%
  summarise(n_individuals = sum(Count), .groups = "drop")

grid <- expand_grid(Species = unique(scol$Species),
                    Lure    = lure_levels,
                    date    = all_dates) %>%
  left_join(distinct(scol, Genus, Species), by = "Species")

summ <- grid %>%
  left_join(summ, by = c("Genus", "Species", "Lure", "date")) %>%
  mutate(n_individuals = replace_na(n_individuals, 0),
         Lure_lab = factor(lure_labels[Lure], levels = lure_labels[lure_levels]))

# species order: abundance-weighted mean collection date (early -> late)
sp_phen <- summ %>%
  group_by(Species) %>%
  summarise(total = sum(n_individuals),
            peak_day = if (sum(n_individuals) > 0)
              as.numeric(weighted.mean(date, n_individuals)) else NA_real_,
            .groups = "drop") %>%
  arrange(peak_day)                # earliest first -> earliest at bottom (y is drawn bottom-up)
summ <- summ %>%
  mutate(Species = factor(Species, levels = sp_phen$Species),
         sp_label = factor(paste0(Species, "  (n=", sp_phen$total[match(Species, sp_phen$Species)], ")"),
                           levels = paste0(sp_phen$Species, "  (n=", sp_phen$total, ")")),
         is_target = Species %in% targets)

write.csv(summ %>%
            select(Genus, Species, Lure, Lure_lab, date, n_individuals) %>%
            arrange(Species, Lure, date),
          file.path(out_data_dir, "scolytinae_species_abundance_by_date_lure.csv"),
          row.names = FALSE)
write.csv(sp_phen, file.path(out_data_dir, "scolytinae_species_phenology_order.csv"),
          row.names = FALSE)

## ---- 4. Heat map: species x date, faceted by lure ------------------------

hm <- summ %>% mutate(
  fill_val = ifelse(n_individuals > 0, log10(n_individuals), NA_real_),
  date_lab = factor(format(date, "%b %d"), levels = format(all_dates, "%b %d")))

brk <- c(1, 3, 10, 30, 100, 300, 1000)
p_hm <- ggplot(hm, aes(date_lab, sp_label, fill = fill_val)) +
  geom_tile(colour = "grey90") +
  geom_text(data = subset(hm, n_individuals >= 10),
            aes(label = n_individuals, colour = n_individuals >= 80), size = 2.4) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  facet_wrap(~ Lure_lab, nrow = 1) +
  scale_fill_viridis_c(option = "mako", direction = -1, na.value = "grey96",
                       breaks = log10(brk), labels = brk,
                       name = "Individuals\n(pooled over sites)") +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Collection date", y = NULL,
       title = "Scolytinae species abundance across the season, by lure",
       subtitle = paste("Rows ordered by abundance-weighted mean date (early-season at bottom).",
                        "Bold = one of the three original target species.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7,
                     face = ifelse(levels(hm$sp_label) %in%
                                   hm$sp_label[hm$is_target], "bold", "plain")),
        strip.text = element_text(face = "bold"),
        legend.key.height = unit(1.2, "cm"))

ggsave(file.path(out_fig_dir, "scolytinae_species_abundance_heatmap.png"),
       p_hm, width = 11, height = 9, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_species_abundance_heatmap.pdf"),
       p_hm, width = 11, height = 9)

## ---- 5. Subfamily-level totals per date x lure --------------------------

tot <- scol %>%
  group_by(Lure, date) %>%
  summarise(n_individuals = sum(Count),
            n_species     = n_distinct(Species), .groups = "drop")
tot <- expand_grid(Lure = lure_levels, date = all_dates) %>%
  left_join(tot, by = c("Lure", "date")) %>%
  mutate(across(c(n_individuals, n_species), ~replace_na(.x, 0)),
         Lure_lab = factor(lure_labels[Lure], levels = lure_labels[lure_levels]))

write.csv(tot %>% select(Lure, Lure_lab, date, n_species, n_individuals),
          file.path(out_data_dir, "scolytinae_subfamily_totals_by_date_lure.csv"),
          row.names = FALSE)

tot_long <- tot %>%
  pivot_longer(c(n_individuals, n_species), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, n_individuals = "Individuals", n_species = "Species richness"),
         metric = factor(metric, levels = c("Individuals", "Species richness")))

p_tot <- ggplot(tot_long, aes(date, value, colour = Lure_lab, group = Lure_lab)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("#4477AA", "#EE6677", "#228833"), name = "Lure") +
  scale_x_date(date_labels = "%b %d", breaks = all_dates) +
  labs(x = "Collection date", y = NULL,
       title = "Scolytinae subfamily: total catch and richness across the season") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_fig_dir, "scolytinae_subfamily_totals_by_lure.png"),
       p_tot, width = 8, height = 6, dpi = 150)
ggsave(file.path(out_fig_dir, "scolytinae_subfamily_totals_by_lure.pdf"),
       p_tot, width = 8, height = 6)

cat("\nDone. Outputs in", out_data_dir, "and", out_fig_dir, "\n")
