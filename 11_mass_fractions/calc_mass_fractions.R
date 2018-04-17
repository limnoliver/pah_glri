# calculate % by mass for each source-sample combination, as well 
# as summaries of all samples



test <- calc_mass_fractions(compound_info = make('samples'),
                            sample_column = 'unique_id',
                            conc_column = 'RESULT',
                            conc_unit = 'ppb',
                            calc_type = 'summary',
                            plot = FALSE)

test2 <- calc_mass_fractions(compound_info = make('samples'),
                            sample_column = 'unique_id',
                            conc_column = 'RESULT',
                            conc_unit = 'ppb',
                            calc_type = 'by_sample',
                            plot = FALSE)

test3 <- calc_mass_fractions(compound_info = make('samples'),
                             sample_column = 'unique_id',
                             conc_column = 'RESULT',
                             conc_unit = 'ppb',
                             calc_type = 'by_sample',
                             plot = TRUE)

quo_sample_column <- sym(sample_column)
quo_conc_column <- sym(conc_column)
quo_compound_column <- sym(compound_column)

library(dplyr)
dat <- filter(compound_info, EPApriority16 %in% TRUE) %>%
  group_by(!!quo_sample_column) %>%
  summarize(sum_EPA16 = sum(!!quo_conc_column)) %>%
  select(quo_sample_column, quo_conc_column)

if (conc_unit == 'ppb') {
  dat <- mutate(dat, conc_mg_kg = ((!!quo_conc_column)/1000))
} else {
  dat <- mutate(dat, conc_mg_kg = !!quo_conc_column)
}


sources <- pah::source_conc %>%
  rename(unique_id = source, concentration = conc_mgkg) %>%
  mutate(type = 'source')

combo <- rbind(samples, sources)

# site-summary analysis
# first take means of sources

source_mean <- group_by(sources, unique_id) %>%
  summarize(concentration = mean(concentration),
            PAHs_used = mean(PAHs_used),
            reference = paste0(reference, collapse = ", ")) %>%
  mutate(units = "mg/kg")

samp_summary <- as.numeric(summary(samples$concentration))

head(source_mean)

mass_fraction <- function(x) {
  temp <- round((x/source_mean$concentration)*100, 2)
  temp <- ifelse(temp <= 100, temp, ">100")
}
source_mean$mass_frac_min <- round((samp_summary[[1]]/source_mean$concentration)*100, 2)
source_mean$mass_frac_Q1 <- round((samp_summary[[2]]/source_mean$concentration)*100, 2)
source_mean$mass_frac_Q2 <- round((samp_summary[[3]]/source_mean$concentration)*100, 2)
source_mean$mass_frac_mean <- round((samp_summary[[4]]/source_mean$concentration)*100, 2)
source_mean$mass_frac_Q3 <- round((samp_summary[[5]]/source_mean$concentration)*100, 2)
source_mean$mass_frac_max <- round((samp_summary[[6]]/source_mean$concentration)*100, 2)

source_mean$mass_frac_min <- ifelse(source_mean$mass_frac_min <= 100, source_mean$mass_frac_min, ">100")
source_mean$mass_frac_Q1 <- ifelse(source_mean$mass_frac_Q1 <= 100, source_mean$mass_frac_Q1, ">100")
source_mean$mass_frac_Q2 <- ifelse(source_mean$mass_frac_Q2 <= 100, source_mean$mass_frac_Q2, ">100")
source_mean$mass_frac_mean <- ifelse(source_mean$mass_frac_mean <= 100, source_mean$mass_frac_mean, ">100")
source_mean$mass_frac_Q3 <- ifelse(source_mean$mass_frac_Q3 <= 100, source_mean$mass_frac_Q3, ">100")
source_mean$mass_frac_max <- ifelse(source_mean$mass_frac_max <= 100, source_mean$mass_frac_max, ">100")

source_mean_wide <- group_by(sources, unique_id) %>%
  summarize(concentration = round(mean(concentration), 2),
            PAHs_used = mean(PAHs_used),
            reference = paste0(reference, collapse = ", ")) %>%
  mutate(units = "mg/kg", type = 'source') %>%
  select(unique_id, concentration)

frac_by_site <- as.data.frame(matrix(ncol = 19, nrow = 1))
names(frac_by_site) <- c("unique_id", source_mean_wide$unique_id)
#frac_by_site[1,1] <- 'source mean concentration'
#frac_by_site[1, 2:19] <- source_mean_wide$concentration

frac_by_site <- bind_rows(frac_by_site, select(samples, unique_id))
frac_by_site$sample_concentration <- c(NA, samples$concentration)

for (i in 2:nrow(frac_by_site)) {
  frac_by_site[i,2:19] <- round((frac_by_site$sample_concentration[i]/frac_by_site[1,2:19])*100, 2)
}

convert_100 <- function(x){
  dat <- ifelse(x <= 100, x, ">100")
  return(dat)
}


#frac_by_site$unique_id <- samples$unique_id
source_names <- source_mean$unique_id
get_possible <- function(x) {
  impossible <- paste0(source_names[which(x > 100)], collapse = ", ")
  unlikely <- paste0(source_names[which(x >= 5 & x <= 100)], collapse = ", ")
  possible <- paste0(source_names[which(x < 5)], collapse = ", ")
   
  return(c(impossible, unlikely, possible))
}
test <- apply(frac_by_site[2:nrow(frac_by_site), 2:19], 1, get_possible)

test <- t(test)
test <- as.data.frame(test)
test$unique_id <- samples$unique_id
names(test)[1:3] <- c('impossible_sources', 'unlikely_sources', 'possible_sources')

frac_by_site_formatted <- as.data.frame(apply(frac_by_site[2:nrow(frac_by_site), 2:19], 2, convert_100))

library(tidyr)
frac_for_plot <- frac_by_site[2:nrow(frac_by_site), 1:19] %>%
  gather(key = "source", value = "mass_fraction", -unique_id)

library(ggplot2)

sample_order <- make('sample_order_pah16')
source_order <- source_mean %>%
  arrange(concentration)

frac_for_plot$source <- factor(frac_for_plot$source, levels = source_order$unique_id)
frac_for_plot$unique_id <- factor(frac_for_plot$unique_id, levels = sample_order)

ggplot(frac_for_plot, aes(x = source, y = unique_id)) +
  geom_tile(aes(fill = cut(mass_fraction, c(-Inf, 5, 100, Inf))), color = 'white') +
  scale_fill_manual(name = 'Mass Fraction', values = c("(-Inf,5]" = '#1a9850',
                                                       "(5,100]" = '#fdae61',
                                                       "(100, Inf]" = '#d73027'), 
                    labels = c('possible (<5%)', 'unlikely (5-100%)', 'impossible (>100%)')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y = "") +
  geom_hline(yintercept = c(18, 36, 54), size = 1)
