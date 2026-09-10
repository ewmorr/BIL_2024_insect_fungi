##############################################################################
# Target-genus phenology -- Ips and Dendroctonus
#
# These EDRR traps are baited specifically for the bark-beetle genera Ips and
# Dendroctonus, so those two lineages are worth looking at on their own rather
# than only as part of the whole-community insect table. This is the first
# script in that line of exploratory work.
#
# Here: number of species and number of individuals caught in each of the two
# target genera, over the six collection dates, faceted by lure. Counts are
# pooled across the four sites (one trap per site x lure), so each point is the
# whole-network catch for that date x lure.
#
# Lure codes in the raw table: E  = Ethanol (generic),
#                              EA = Alpha-pinene + Ethanol (pine bark beetles),
#                              I  = Ips lure (ipsenol / ipsdienol).
#
# Data source: data/2024_insect_data/EDRRPatho_NHTraps_InsectData_MFD.csv
# (raw long-format trap catch -- same file insect_dat_process.R starts from).
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

out_fig_dir  <- "figures/insect_exploratory/target_genera"
out_data_dir <- "data/insect_exploratory/target_genera"
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

target_genera <- c("Ips", "Dendroctonus")

lure_levels <- c("E", "EA", "I")
lure_labels <- c(E  = "Ethanol",
                 EA = "Alpha-pinene + Ethanol",
                 I  = "Ips lure")

## ---- 1. Load + tidy the raw long-format trap catch --------------------------

dat <- read.csv("data/2024_insect_data/EDRRPatho_NHTraps_InsectData_MFD.csv",
                check.names = FALSE)
names(dat) <- make.names(names(dat))

# Trailing all-blank rows in the sheet.
dat <- dat %>% filter(!is.na(Trap.Code), Trap.Code != "")
dat[dat == ""] <- NA

dat <- dat %>%
  mutate(date = mdy(Collection.Date),
         Count = suppressWarnings(as.numeric(Count)))

# One record (P / Ips lure / 6/12/24) is blank apart from its trap code -- a
# serviced trap with no catch entered. It carries no genus, so dropping it does
# not affect the target-genus counts.
n_blank <- sum(is.na(dat$Count))
if (n_blank > 0) {
  cat("Dropping", n_blank, "catch record(s) with no Count entered (no genus attached).\n")
  dat <- dat %>% filter(!is.na(Count))
}
stopifnot(!any(is.na(dat$date)))

## ---- 2. Restrict to the two target genera ----------------------------------

tg <- dat %>% filter(Genus %in% target_genera)

cat("Target-genus records:\n")
tg %>% count(Genus, Species, wt = Count, name = "individuals") %>% arrange(Genus, desc(individuals)) %>% print()
cat("\nGenus x lure (records):\n"); print(table(tg$Genus, tg$Lure))

# Every target-genus row here is resolved to species (no genus-only "Ips sp."
# records in this dataset), so "number of species" = distinct Species values.
stopifnot(all(!is.na(tg$Species)))

## ---- 3. Pool across sites: per date x lure x genus -------------------------

summ <- tg %>%
  group_by(Genus, Lure, date) %>%
  summarise(n_individuals = sum(Count),
            n_species     = n_distinct(Species),
            .groups = "drop")

# Fill the date x lure x genus grid so lines drop to zero on dates with no catch.
grid <- expand_grid(Genus = target_genera,
                    Lure  = lure_levels,
                    date  = sort(unique(dat$date)))

summ <- grid %>%
  left_join(summ, by = c("Genus", "Lure", "date")) %>%
  mutate(across(c(n_individuals, n_species), ~replace_na(.x, 0)),
         Genus     = factor(Genus, levels = target_genera),
         Lure_lab  = factor(lure_labels[Lure], levels = lure_labels[lure_levels]))

write.csv(summ %>% select(Genus, Lure, Lure_lab, date, n_species, n_individuals),
          file.path(out_data_dir, "target_genus_catch_by_date_lure.csv"),
          row.names = FALSE)

## ---- 4. Plot: metric (rows) x lure (cols), coloured by genus --------------

plot_df <- summ %>%
  pivot_longer(c(n_species, n_individuals),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         n_species     = "Number of species",
                         n_individuals = "Number of individuals"),
         metric = factor(metric, levels = c("Number of individuals", "Number of species")))

genus_cols <- c(Ips = "#1b9e77", Dendroctonus = "#d95f02")

p <- ggplot(plot_df, aes(date, value, colour = Genus, group = Genus)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_grid(metric ~ Lure_lab, scales = "free_y", switch = "y") +
  scale_colour_manual(values = genus_cols, name = "Genus") +
  scale_x_date(date_labels = "%b %d", breaks = sort(unique(dat$date))) +
  labs(x = "Collection date", y = NULL,
       title = "Target bark-beetle genera: catch over the season by lure",
       subtitle = "Pooled across the four sites (one trap per site x lure)") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "grey92", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_fig_dir, "target_genus_phenology_by_lure.png"),
       p, width = 9, height = 6, dpi = 150)
ggsave(file.path(out_fig_dir, "target_genus_phenology_by_lure.pdf"),
       p, width = 9, height = 6)

## ---- 5. Species x lure breakdown -----------------------------------------
# Only three target species are present (Ips grandicollis, Ips pini,
# Dendroctonus valens), so "number of species" is no longer a useful axis
# here -- this section is individuals per species, over date, species x lure.

sp_order <- tg %>% count(Genus, Species, wt = Count, name = "n") %>%
  arrange(Genus, desc(n)) %>% pull(Species)

summ_sp <- tg %>%
  group_by(Genus, Species, Lure, date) %>%
  summarise(n_individuals = sum(Count), .groups = "drop")

grid_sp <- expand_grid(Species = sp_order,
                       Lure    = lure_levels,
                       date    = sort(unique(dat$date))) %>%
  left_join(distinct(tg, Genus, Species), by = "Species")

summ_sp <- grid_sp %>%
  left_join(summ_sp, by = c("Genus", "Species", "Lure", "date")) %>%
  mutate(n_individuals = replace_na(n_individuals, 0),
         Species  = factor(Species, levels = sp_order),
         Lure_lab = factor(lure_labels[Lure], levels = lure_labels[lure_levels]))

write.csv(summ_sp %>% select(Genus, Species, Lure, Lure_lab, date, n_individuals) %>%
            arrange(Species, Lure, date),
          file.path(out_data_dir, "target_species_catch_by_date_lure.csv"),
          row.names = FALSE)

# Season totals per species x lure (the "how well does each lure catch each
# species" summary).
tot_sp <- summ_sp %>%
  group_by(Genus, Species, Lure, Lure_lab) %>%
  summarise(n_individuals = sum(n_individuals), .groups = "drop")
write.csv(tot_sp %>% select(Genus, Species, Lure, Lure_lab, n_individuals) %>%
            arrange(Species, Lure),
          file.path(out_data_dir, "target_species_totals_by_lure.csv"),
          row.names = FALSE)
cat("\nSeason totals, individuals per species x lure:\n")
print(tot_sp %>% select(Species, Lure, n_individuals) %>%
        pivot_wider(names_from = Lure, values_from = n_individuals))

sp_cols <- c("Ips grandicollis"    = "#1b9e77",
             "Ips pini"            = "#66c2a5",
             "Dendroctonus valens" = "#d95f02")

# 5a. Phenology: individuals over date, species (rows) x lure (cols)
p_sp <- ggplot(summ_sp, aes(date, n_individuals, colour = Species, group = Species)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_grid(Species ~ Lure_lab, scales = "free_y", switch = "y") +
  scale_colour_manual(values = sp_cols, guide = "none") +
  scale_x_date(date_labels = "%b %d", breaks = sort(unique(dat$date))) +
  labs(x = "Collection date", y = "Number of individuals",
       title = "Target bark-beetle species: catch over the season by lure",
       subtitle = "Pooled across the four sites; note the free y-axis per species") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left = element_text(angle = 0, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(out_fig_dir, "target_species_phenology_by_lure.png"),
       p_sp, width = 9, height = 6, dpi = 150)
ggsave(file.path(out_fig_dir, "target_species_phenology_by_lure.pdf"),
       p_sp, width = 9, height = 6)

# 5b. Season totals: individuals per species x lure
p_tot <- ggplot(tot_sp, aes(Lure_lab, n_individuals, fill = Species)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_individuals), vjust = -0.3, size = 3) +
  facet_wrap(~ Species, nrow = 1) +
  scale_fill_manual(values = sp_cols, guide = "none") +
  labs(x = "Lure", y = "Number of individuals (whole season, all sites)",
       title = "Target bark-beetle species: total catch by lure") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(out_fig_dir, "target_species_totals_by_lure.png"),
       p_tot, width = 8, height = 4.2, dpi = 150)
ggsave(file.path(out_fig_dir, "target_species_totals_by_lure.pdf"),
       p_tot, width = 8, height = 4.2)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
