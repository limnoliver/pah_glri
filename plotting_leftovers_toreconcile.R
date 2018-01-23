


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






