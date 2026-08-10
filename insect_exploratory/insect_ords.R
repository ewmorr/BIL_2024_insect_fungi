library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
source("library/library.R")

sp_tab = read.csv("data/2024_insect_data/insect_species_tab.csv")
head(sp_tab)
str(sp_tab)
# there is a missing value with all zeros 
sp_tab = sp_tab[!is.na(sp_tab$Finest.ID),]

# set rownames to finest.ID for t() 
rownames(sp_tab) = sp_tab$Finest.ID

sp_tab

sp_tab.t = t(sp_tab %>% select(where(is.numeric)))
head(sp_tab.t)

rowSums(sp_tab.t)
sum(rowSums(sp_tab.t))

sp_nmds = metaMDS(sp_tab.t[rowSums(sp_tab.t) > 0,], distance = "bray", binary = F, try = 20, trymax = 100)
metaMDS(sp_tab.t[rowSums(sp_tab.t) > 0,], distance = "bray", binary = F, try = 20, trymax = 100, previous.best = sp_nmds)
# not getting a repeatable solution. Let's look at colSums
plot(sort(colSums(sp_tab.t))) 
# lot's of things not nrepped well
plot(sort(colSums(sp_tab %>% filter(Order == "Coleoptera") %>% select(where(is.numeric)) %>% t() ))) 

sp_tab.t[sp_tab.t > 0] = 1

plot(sort(colSums(sp_tab.t)))

# transformations don't help
metaMDS(sqrt(sp_tab.t[rowSums(sp_tab.t) > 0,]), distance = "bray", binary = F, try = 20, trymax = 100)
metaMDS(log10(sp_tab.t[rowSums(sp_tab.t) > 0,]+1), distance = "bray", binary = F, try = 20, trymax = 100)

# Coleoptera is by far most abundant order (~100x greater than next most)
metaMDS(sp_tab.coleop_genus.t[sp_tab.coleop_genus.t], distance = "bray", binary = F, try = 20, trymax = 100)
#still no love

##########################################################
# lot's of things occur once. We have a couple of options. 
# First sum at a higher tax level
# second filter low occurrence

# we already did some looks at the abundance of different families
# we probably want to filter to scolytid since these traps are targeting
# them and they are most abundant (or Curculionidae sub family)
##############################

# Genus
sp_tab.coleop_genus = sp_tab %>% 
    filter(Order == "Coleoptera") %>%
    summarise(.by = c(Class, Order, Family, Subfamily, Genus), across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

sp_tab.coleop_genus.t = t(sp_tab.coleop_genus %>% select(where(is.numeric)))
plot(sort(colSums(sp_tab.coleop_genus.t)))
sp_tab.coleop_genus.t[sp_tab.coleop_genus.t > 0] = 1
plot(sort(colSums(sp_tab.coleop_genus.t)))

# Family
sp_tab.coleop_genus = sp_tab %>% 
    filter(Order == "Coleoptera") %>%
    summarise(.by = c(Class, Order, Family), across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

sp_tab.coleop_genus.t = t(sp_tab.coleop_genus %>% select(where(is.numeric)))
plot(sort(colSums(sp_tab.coleop_genus.t)))
sp_tab.coleop_genus.t[sp_tab.coleop_genus.t > 0] = 1
plot(sort(colSums(sp_tab.coleop_genus.t)))

rowSums(sp_tab.coleop_genus.t) == 0
# only 8 families only occur once


# The Latridiidae (1872) and Curculionidae (5729) are by far the most abundant fams in Coleoptera
# (rest are under 400)

sp_tab.sums = data.frame(
    sp_tab %>%
        select(!where(is.numeric)),
    total.inds = sp_tab %>%
        select(where(is.numeric)) %>%
        rowSums()
)

family_sums = sp_tab.sums %>%
    summarize(.by = c(Class, Order, Family), across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    arrange(total.inds)

total_inds = sum(family_sums$total.inds)
#10926 individuals
family_sums %>% 
    filter(Family %in% c("Curculionidae", "Latridiidae")) %>%
    select(total.inds) %>% sum()
7601/10926
# 69% of total fall into these families. Nice.

#####################################################
#####################################################
# try filtering to just Latridiidae and Curculionidae 
sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t = t(sp_tab.curcus_latris %>% select(where(is.numeric)))
length(sp_tab.curcus_latris$Finest.ID) == length(sp_tab.curcus_latris$Finest.ID %>% unique())
colnames(sp_tab.curcus_latris.t) = sp_tab.curcus_latris$Finest.ID

colSums(sp_tab.curcus_latris.t)
# remove singletons
#52 of 63 remaining
sp_tab.curcus_latris.t[,colSums(sp_tab.curcus_latris.t) > 1] -> sp_tab.curcus_latris.no_singleton
rowSums(sp_tab.curcus_latris.no_singleton) # there is a sample with no inds
sp_tab.curcus_latris.no_singleton[rowSums(sp_tab.curcus_latris.no_singleton) > 0,] -> sp_tab.curcus_latris.no_singleton

# now try NMDS
metaMDS(sp_tab.curcus_latris.no_singleton, distance = "bray", binary = F, try = 20, trymax = 100)
# still high stress but we get a high repetition
metaMDS(sp_tab.curcus_latris.no_singleton, distance = "bray", binary = T, try = 20, trymax = 100)

# oh we get low stress without the transform
# best solution repeated once in 20 tries. acceptable
metaMDS(sp_tab.curcus_latris.no_singleton, distance = "bray", binary = F, try = 20, trymax = 100, autotransform = F)
# let's try the full dataset with no autotrans
metaMDS(sp_tab.t[rowSums(sp_tab.t) > 0,], distance = "bray", binary = F, try = 20, trymax = 100, autotransform = F)
# stress again at >0.23 and best only repeats once in 100

metaMDS(log10(sp_tab.curcus_latris.no_singleton+1), distance = "bray", binary = F, try = 20, trymax = 100, autotransform = F)
metaMDS(sqrt(sp_tab.curcus_latris.no_singleton), distance = "bray", binary = F, try = 20, trymax = 100, autotransform = F)
# funny, the transforms increases stress. Probably because we are downweighting the dominants
# let's roll with the untransformed data

#################################################
# run the Curculionidae and Latridiidae only NMDS
curcus_latris.nmds = metaMDS(sp_tab.curcus_latris.no_singleton, distance = "bray", binary = F, try = 20, trymax = 100, autotransform = F)
plot(curcus_latris.nmds)

saveRDS(curcus_latris.nmds, "data/2024_insect_data/Curculionidae_and_Latridiidae.NMDS.RDS")

# join site data with metadata
curcus_latris.nmds.sites = data.frame(scores(curcus_latris.nmds)$sites)
curcus_latris.nmds.species = data.frame(scores(curcus_latris.nmds)$species)

# read the sequence metadata and filter down to insect samples on fisrt read through
#metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.redo.csv")
#row.names(sp_tab.curcus_latris.no_singleton) %>% sub("\\.\\.", " ", .) %>% gsub("\\.", "-", .) -> table_ids
#metadata %>% filter(sampleID %in% table_ids) %>% nrow()
#table_ids[!table_ids %in% metadata$sampleID]
#metadata = metadata %>% filter(sampleID %in% table_ids)
#nrow(curcus_latris.nmds.sites)
#metadata %>% summarize(.by = sampleID, n = n())
#write.csv(metadata, "data/2024_metadata/insect_community_metadata.csv", row.names = F)
metadata = read.csv("data/metadata/insect_community_metadata.csv")

curcus_latris.nmds.sites$sampleID = row.names(curcus_latris.nmds.sites) %>% sub("\\.\\.", " ", .) %>% gsub("\\.", "-", .) 
curcus_latris.nmds.sites.metadata = left_join(curcus_latris.nmds.sites, metadata)

curcus_latris.nmds.sites.metadata$mdy = lubridate::mdy(curcus_latris.nmds.sites.metadata$CollectionDate)

# run PERMANOVA
mod1 = adonis2(sp_tab.curcus_latris.no_singleton ~ Risk + Risk:Site + Risk:Site:Lure, data = curcus_latris.nmds.sites.metadata, by = "terms")
mod1

# To correctly calculate p-value for an effect of site risk we need to permute
# the whole site-level block of samples and not shuffle between sites.
# However, we also only have two sites per each of the risk categories
# so this is probably extraneous for this dataset (we can calc for the fungal
# communities for point of compairison)
# Let's plan to look at strength of effects across the entire fungal community data set
# and within these sites only
# and also calculate a measure of concordance between insect/fungal comms)
# 
# We need to find out if we can compare to the spore trap dataset directly
# or if there is a way to calc concordance bweteen our data and that

# Builiding correct calls for adoni2 here
# # date
perm_date <- how(
  blocks = curcus_latris.nmds.sites.metadata$Site,
  plots  = Plots(strata = curcus_latris.nmds.sites.metadata$trapID, type = "none"),
  within = Within(type = "free"),
  nperm  = 999
)

adonis2(sp_tab.curcus_latris.no_singleton ~ Site + Lure + mdy, data = curcus_latris.nmds.sites.metadata, permutations = perm_date, by = "margin")

#         Df SumOfSqs      R2       F Pr(>F)    
#Site      3   2.0162 0.09739  3.3367  0.001 ***
#Lure      2   2.9315 0.14160  7.2771  0.001 ***
#mdy       1   2.9332 0.14168 14.5627  0.001 ***
#Residual 64  12.8908 0.62265                   
#Total    70  20.7032 1.00000                   

adonis2(sp_tab.curcus_latris.no_singleton ~ Site/Lure + Site/Lure/mdy, data = curcus_latris.nmds.sites.metadata, permutations = perm_date, by = "margin")


# Lure

# we have unequal samples for each trap bc a date is missing so we subsample
set.seed(1)

trap_counts <- table(curcus_latris.nmds.sites.metadata$trapID)
min_n <- min(trap_counts)

keep_rows <- unlist(lapply(split(seq_len(nrow(curcus_latris.nmds.sites.metadata)), curcus_latris.nmds.sites.metadata$trapID), function(idx) {
  sample(idx, size = min_n)
}))

meta_bal      <- curcus_latris.nmds.sites.metadata[keep_rows, ]
sp_matrix_bal <- sp_tab.curcus_latris.no_singleton[keep_rows, , drop = FALSE]

table(meta_bal$trapID)  # confirm equal now
table(meta_bal$Site)  # confirm equal now

perm_lure <- how(
  blocks = meta_bal$Site,
  plots  = Plots(strata = meta_bal$trapID, type = "free"),
  within = Within(type = "none"),
  nperm  = 999
)


adonis2(sp_matrix_bal ~ Site + Lure + mdy, data = meta_bal, permutations = perm_lure, by = "margin")

#         Df SumOfSqs      R2      F Pr(>F)    
#Site      3   1.6869 0.09938 2.3898  0.001 ***
#Lure      2   2.5819 0.15210 5.4866  0.002 ** 
#Residual 54  12.7059 0.74852                  
#Total    59  16.9748 1.00000               

# the lure r2 here is affected by removing random samples to balance
# use the lure r2 from the full model but use this p-val for lure

# strong and almost equal effects of Site and Date
# not surprising but good to find

#adonis2(formula = sp_tab.curcus_latris.no_singleton ~ Site + Site:Lure + Site:mdy, data = curcus_latris.nmds.sites.metadata, by = "terms")
#          Df SumOfSqs      R2      F Pr(>F)    
#Site       3   2.0590 0.09945 3.6245  0.001 ***
#Site:Lure  8   4.1608 0.20098 2.7466  0.001 ***
#Site:mdy   4   4.0682 0.19650 5.3709  0.001 ***
#Residual  55  10.4151 0.50307                  
#Total     70  20.7032 1.00000                  
#---

ggplot(curcus_latris.nmds.sites.metadata, aes(x = NMDS1, y = NMDS2, fill = Lure, shape = Site)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(21,22,23,24)) +
    scale_fill_brewer(palette = "Dark2") +
    guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2, override.aes = list(shape = 22))) +
    annotate(geom = "text", label = "stress = 0.14", x = 0.65, y = -2.1) +
    theme_bw() -> p1
p1
#date_breaks(as.integer(curcus_latris.nmds.sites.metadata %>% filter(!is.na(mdy)) %>% pull(mdy)))
#five_cols_gradient_palette = c('#ca0020','#ef8a62','#f7f7f7','#67a9cf','#2166ac')

ggplot(curcus_latris.nmds.sites.metadata, aes(x = NMDS1, y = NMDS2, fill = mdy, shape = Site)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(21,22,23,24)) +
    annotate(geom = "text", label = "stress = 0.14", x = 0.65, y = -2.1) +
    theme_bw() -> p2
p2

pdf("figures/FEDRR_all_2024/insects.date_and_lure.NH.pdf", width = 9, height = 7)
p1
p2
dev.off()



#################################################
# Also testing an NMDS with everything but these


head(metadata)




head(metadata)
nrow(metadata)
nrow(sp_tab.curcus_latris.no_singleton)
length(table_ids)
table_ids[!table_ids %in% metadata$sampleID]
