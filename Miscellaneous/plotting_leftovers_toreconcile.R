# plot source chi squareds across all samples
# first show dropping creosote

# order sources by median value
creosote <- unique(grep('creosote', all.summary$source, ignore.case = T, value = T))
all.no.creosote <- filter(all.summary, !(source %in% creosote))
order.vals <- all.no.creosote %>%
  group_by(source) %>%
  summarize(med = median(sum_chi2)) %>%
  arrange(med)

all.no.creosote$source <- factor(all.no.creosote$source, levels = order.vals$source)

p <- ggplot(all.no.creosote, aes(x = source, y = sum_chi2)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = '', y = 'Sum Chi2')

ggsave('6_profile_data/doc/chi2_allsites_nocreosote.png', p)

all.no.creosote <- left_join(all.no.creosote, samples_16)

p <- p +
  facet_grid(~Priority16_bin)

ggsave('6_profile_data/doc/chi2_allsites_binnedbypriority16.png', p)

## calculate what is source with the smallest distance by each site
which.source <- group_by(all.no.creosote, sample_id) %>%
  summarize(min_source = source[which.min(sum_chi2)])

compound.drop <- unique(all.profs$casrn[is.na(all.profs$source_prop_conc)])
creosote.summary <- filter(all.profs, casrn != compound.drop) %>%
  group_by(sample_id, source) %>%
  summarize(sum_chi2 = sum(chi2))

order.vals <- creosote.summary %>%
  group_by(source) %>%
  summarize(med = median(sum_chi2)) %>%
  arrange(med)
creosote.summary$source <- factor(creosote.summary$source, levels = order.vals$source)

p <- ggplot(creosote.summary, aes(x = source, y = sum_chi2)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = '', y = 'Sum Chi2')

ggsave('6_profile_data/doc/chi2_allsites_withcreosote.png', p)

# create plot of proportional concentrations across all samples compared to individual sources
dat <- out[[1]]
sample.mean <- group_by(dat, Compound) %>%
  summarize(mean_prop_conc = mean(prop_conc),
            sd_prop_conc = sd(prop_conc)) %>%
  rename(sample_prop_conc = mean_prop_conc)

sources <- select(dat, Compound, source, source_prop_conc) %>%
  rename(prop_conc = source_prop_conc, profile = source) %>%
  distinct() %>%
  filter(!(profile %in% c('Creosote_railway_ties', 'Creosote_product')))

profiles <- left_join(sources, sample.mean) %>%
  left_join(distinct(dat[,c('Compound', 'molwt')])) %>%
  filter(profile %in% c('CT_dust_7', 'Vehicle_traffic_avg', 'Tire_particles', 
                        'Residential_heating', 'Asphalt', 'Pine_combustion_1'))

#################
# plots from double ratios

# calculate euclidean distances
library(tidyr)
library(dplyr)


# samp1 <- select(ratios, sample_id, sample_type, Anth_AnthPhen) %>%
#   filter(sample_type == 'sample') %>%
#   select(-sample_type) %>%
#   rename(Anth_AnthPhen_sample = Anth_AnthPhen)

structure1 <- select(ratios, sample_id, sample_type) %>%
  filter(sample_type == 'sample') %>%
  select(sample_id)

structure2 <- select(ratios, sample_id, sample_type) %>%
  filter(sample_type == 'source') %>%
  select(sample_id)

structure <- expand.grid(structure1[[1]], structure2[[1]])

ratios_sub <- select(ratios, -sample_type)

comp <- rename(structure, sample = Var1, source = Var2) %>%
  left_join(ratios_sub, by = c('sample' = 'sample_id')) %>%
  rename(Anth_AnthPhen_sam = Anth_AnthPhen,
         Flua_FluaPyr_sam = Flua_FluaPyr,
         Baa_BaaCh_sam = Baa_BaaCh,
         Indpy_IndpyBghip_sam = Indpy_IndpyBghip) %>%
  left_join(ratios_sub, by = c('source' = 'sample_id')) %>%
  rename(Anth_AnthPhen_src = Anth_AnthPhen,
         Flua_FluaPyr_src = Flua_FluaPyr,
         Baa_BaaCh_src = Baa_BaaCh,
         Indpy_IndpyBghip_src = Indpy_IndpyBghip)

comp <- comp %>%
  mutate(diff1 = sqrt((Anth_AnthPhen_sam - Anth_AnthPhen_src)^2 +
                        (Flua_FluaPyr_sam - Flua_FluaPyr_src)^2),
         diff2 = sqrt((Indpy_IndpyBghip_sam - Indpy_IndpyBghip_src)^2 +
                        (Flua_FluaPyr_sam - Flua_FluaPyr_src)^2),
         diff3 = sqrt((Indpy_IndpyBghip_sam - Indpy_IndpyBghip_src)^2 + 
                        (Baa_BaaCh_sam - Baa_BaaCh_src)^2))
ranks <- comp %>%
  group_by(sample) %>%
  mutate(rank1 = row_number(diff1),
         rank2 = row_number(diff2), 
         rank3 = row_number(diff3)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(sumranks = sum(rank1, rank2, rank3, na.rm = T))

# find how many times each source is the top source by sample

by.sources.rank.counts <- ranks %>%
  group_by(source) %>%
  summarize(top_diff1 = ifelse(anyNA(rank1), NA, length(which(rank1 %in% 1))),
            top_diff2 = ifelse(anyNA(rank2), NA, length(which(rank2 %in% 1))),
            top_diff3 = ifelse(anyNA(rank3), NA, length(which(rank3 %in% 1))))

by.sources.rank.counts$n_poss <- (3-rowSums(is.na(by.sources.rank.counts)))*71
by.sources.rank.counts$n_prop <- (rowSums(by.sources.rank.counts[,c('top_diff1', 'top_diff2', 'top_diff3')], na.rm = T)/by.sources.rank.counts$n_poss)*100
by.sources.rank.counts <- filter(by.sources.rank.counts, !is.na(n_prop))
by.sources.rank.counts <- mutate(by.sources.rank.counts, 
                                 percent_top = ifelse(n_prop == 0, 0.1, n_prop))

mean.sources <- ranks %>%
  select(source, diff1, diff2, diff3, rank1, rank2, rank3) %>%
  group_by(source) %>%
  summarize_all(mean) %>%
  ungroup() %>%
  mutate_at(vars(diff1, diff2, diff3), funs(round(., digits = 2))) %>%
  mutate_at(vars(rank1, rank2, rank3), funs(round(., digits = 0))) %>%
  rename(FFP_AAP_meandiff = diff1, FFP_IIB_meandiff = diff2, BBC_IIB_meandiff = diff3,
         FFP_AAP_meanrank = rank1, FFP_IIB_meanrank = rank2, BBC_IIB_meanrank = rank3)

mean.sources <- left_join(mean.sources, by.sources.rank.counts, by = 'source') 
mean.sources <- mean.sources %>%
  rename(FFP_AAP_n_toprank = top_diff1, FFP_IIB_n_toprank = top_diff2, BBC_IIB_n_toprank = top_diff3)

################# plot means by sources
# this should only be done in plotting scenarios, not main data processing
long_diff <- select(mean.sources, source:BBC_IIB_meandiff) %>%
  gather(double_ratio, double_ratio_meandiff, -source) %>%
  mutate(double_ratio = gsub('_meandiff', '', double_ratio))

long_rank <- select(mean.sources, source, FFP_AAP_meanrank:BBC_IIB_meanrank) %>%
  gather(double_ratio, double_ratio_meanrank, -source) %>%
  mutate(double_ratio = gsub('_meanrank', '', double_ratio))

mean.sources.long <- left_join(long_diff, long_rank, by = c('source', 'double_ratio'))

order.sources <- filter(mean.sources.long, double_ratio == "FFP_AAP") %>%
  arrange(double_ratio_meandiff)

mean.sources.long$source <- as.factor(mean.sources.long$source)
mean.sources.long$source <- factor(mean.sources.long$source, levels = order.sources[[1]])
mean.sources.long$double_ratio <- as.factor(mean.sources.long$double_ratio)
levels(mean.sources.long$double_ratio) <- c('BaA/(BaA+Ch) : IndPy/(IndPy+BghiP)',
                                            'FluA/(FluA+Pyr) : Anth/(Anth+Phen)',
                                            'FluA/(FluA+Pyr) : IndPy/(IndPy+BghiP)')
mean.sources.long$double_ratio <- factor(mean.sources.long$double_ratio, levels(mean.sources.long$double_ratio)[c(2,3,1)])

p <- ggplot(mean.sources.long, aes(x = source, y = double_ratio_meandiff)) +
  geom_bar(stat = 'identity', position = 'dodge', aes(fill = double_ratio_meanrank), color = 'black', size = 0.4) +
  scale_fill_gradient2(low = '#67001f', mid = '#f7f7f7', high = '#053061', 
                       midpoint = median(mean.sources.long$double_ratio_meanrank, na.rm = T),
                       name = "Mean Rank") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), strip.background = element_blank()) +
  facet_wrap(~double_ratio, nrow = 3) +
  labs(x = 'Source', y = 'Mean distance between source and samples')

ggsave('7_double_ratios/doc/mean_distance_barchart_withrank.png', p, height = 5, width = 8)

####################################################
p <- ggplot(by.sources.rank.counts, aes(x = reorder(source, n_prop), y = n_prop)) +
  geom_bar(stat = 'identity', aes(fill = factor(n_poss))) +
  scale_fill_brewer(palette = 'Dark2', type = 'qual', name = 'Number of Comparisons') +
  labs(x = 'Sources', y = 'Percent Times Top Source in\nDouble Ratio Comparisons') +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = c(0,1), legend.justification = c(0,1),
        legend.box.background = element_rect(color = 'black'), 
        legend.box.margin = margin(rep(1,4))) +
  guides(fill = guide_legend(nrow = 1))

ggsave('7_double_ratios/doc/percent_top_source.png', height = 4, width = 7)  

############# by samples
by.samples <- ranks %>%
  group_by(sample) %>%
  summarize(top_source1 = source[which.min(diff1)],
            top_source2 = source[which.min(diff2)],
            top_source3 = source[which.min(diff3)])

top_sources <- by.sources$source[by.sources$n_prop > 5]
top_sources <- na.omit(top_sources)

##################################
# plot top sources by sample
# use this code only for plotting
by.samples.long <- by.samples %>%
  gather(key = comparison, value = top_var, -sample) %>%
  mutate(top_var_cat = ifelse(top_var %in% top_sources, top_var, 'other')) %>%
  left_join(prepped_totals[,c('sample_id', 'Priority16')], by = c('sample' = 'sample_id'))

p <- ggplot(by.samples.long, aes(x = comparison, y = reorder(sample, Priority16))) +
  geom_tile(aes(fill = top_var_cat), color = 'white') +
  scale_fill_manual(values = c('#8c510a', '#bf812d', '#c51b7d',
                               '#01665e', '#35978f', 'gray', '#b2182b'), 
                    name = "Source") +
  labs(x = "", y = "")

head(by.samples.long)


ggsave('double_ratio_top_sources.png', p, height = 8, width = 5)
# calculate overall min distance
# this would mean grouping by source, and finding the mean distance for
# each plot or summing three plots ED and finding min(mean)
temp<-sqrt((samplesC$IP_IPBghiP[i]-sourcesC$IP_IP.Bghi)^2+(samplesC$BaA_BaACh[i]-sourcesC$BaA_228)^2)

# rank each sample for all sources,
# sum the ranks to find the closest source-sample combination





