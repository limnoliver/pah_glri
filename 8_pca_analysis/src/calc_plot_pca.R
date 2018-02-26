sample_profiles <- make('profiles')
sample_profiles <- sample_profiles$profiles
source_profiles <- pah::source_profiles
library(dplyr)
library(tidyr)
library(stats)
library(ggplot2)

# get source profiles transposed and ready for pca
source_profiles_t <- select(source_profiles, -Compound, -pcode, -casrn, -molwt, -Abbreviation) %>%
  t()

source_profiles_t <- data.frame(source_profiles_t)
names(source_profiles_t) <- source_profiles$Abbreviation
source_profiles_t$type <- 'source'

# get sample profiles transposed and ready for pca
sample_profiles_t <- select(sample_profiles, sample_id, Abbreviation, prop_conc) %>%
  distinct() %>%
  spread(key = Abbreviation, value = prop_conc)

sample_profiles_t <- as.data.frame(sample_profiles_t)
row.names(sample_profiles_t) <- sample_profiles_t$sample_id
sample_profiles_t <- select(sample_profiles_t, -sample_id)

sample_profiles_t$type <- 'sample'
# put together
profiles_t <- rbind(source_profiles_t, sample_profiles_t)
profiles_t$sample_id <- row.names(profiles_t)

# drop any rows with NA values
profiles_t <- drop_na(profiles_t)

# drop any rows (samples) with 0 values (BDL)
profiles_t <- filter_all(profiles_t, all_vars(. != 0))

sample_id <- profiles_t$sample_id
type <- profiles_t$type

profiles_t <- select(profiles_t, -type, -sample_id)
profiles_t <- data.frame(profiles_t)
row.names(profiles_t) <- sample_id

# run pca
pca <- prcomp(profiles_t, center = T, scale = T)
pca_out <- summary(pca)$importance

# determine what dimensions the plotting space should be in
n_pca <- length(which(pca_out[2,] >= 0.10))
if (n_pca == 1) {
  n_pca <- 2
} else {
  n_plots <- factorial(n_pca)/(factorial(2)*factorial(n_pca - 2))
}


# get data in format to be used with ggplot
pca_plot_dat <- data.frame(pca$x[,1:n_pca])
pca_plot_dat$sample_id <- row.names(pca_plot_dat)
pca_plot_dat$type <- type

# plot
library(GGally)
library(ggrepel)
library(gridExtra)

plot_num <- 1
plot_list <- list()
my_labels <- pca_plot_dat$sample_id
my_labels[pca_plot_dat$type == 'sample'] <- ""

for(xcol in 1:3){
  for(ycol in 2:4){
    if (xcol >= ycol) {next}
    
    temp_dat <- data.frame(x = pca_plot_dat[,xcol],
                           y = pca_plot_dat[,ycol], 
                           type = pca_plot_dat$type)
   
    p <- ggplot(data = temp_dat, aes(x = x, y = y)) +
        geom_point(aes_string(color = 'type'), alpha = 0.5, show.legend = F) +
        geom_text_repel(data = temp_dat,
              aes(x = x, y = y, label = my_labels)) +
        scale_color_manual(values = c('black', 'red')) +
        labs(x = paste('Component', xcol), y = paste('Component', ycol)) +
        theme_classic()
    
    plot_list[[plot_num]] <- p
    plot_num <- plot_num + 1
  }
}

if (n_pca < 4) {
  p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
} else if(n_pca == 4) {
  p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 3)
} else if (n_pca > 4) {
  p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 5)
}

# now calculate euclidean distances and create boxplot
structure1 <- select(pca_plot_dat, sample_id, type) %>%
  filter(type == 'sample') %>%
  select(sample_id)

structure2 <- select(pca_plot_dat, sample_id, type) %>%
  filter(type == 'source') %>%
  select(sample_id)

structure <- expand.grid(structure1[[1]], structure2[[1]])

pca_plot_dat_sub <- select(pca_plot_dat, -type)

# left join by sample ids and rename columns
comp <- rename(structure, sample = Var1, source = Var2) %>%
  left_join(pca_plot_dat_sub, by = c('sample' = 'sample_id'))

sample_colnames <- paste0('PC', 1:n_pca, '_sam')
names(comp)[3:(2+n_pca)] <- sample_colnames

#left join by source ids and rename colums
comp <- comp %>%
  left_join(pca_plot_dat_sub, by = c('source' = 'sample_id'))

source_colnames <- paste0('PC', 1:n_pca, '_src')
names(comp)[(3+n_pca):ncol(comp)] <- source_colnames

# find squared differences between source and sample of 
# each component to reduce dimensions
pca_names <- names(pca_plot_dat)[1:n_pca]
comp_diff <- data.frame(sample = comp$sample, 
                        source = comp$source)
for (i in 1:n_pca) {
  temp_col <- grep(pca_names[i], names(comp))
  temp_diff <- (comp[,temp_col[1]] - comp[,temp_col[2]])^2
  comp_diff[,2+i] <- temp_diff
}

comp_diff <- comp_diff %>%
  mutate(euc_dist = sqrt(rowSums(.[3:(2+n_pca)]))) %>%
  select(sample, source, euc_dist)

# boxplot
median_diff <- comp_diff %>%
  group_by(source) %>%
  summarize(median = median(euc_dist)) %>%
  arrange(median)

comp_diff$source <- factor(comp_diff$source, levels = median_diff$source)


ggplot(comp_diff, aes(x = source, y = euc_dist)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Euclidean distance\n(zero = identical to sample)", x = "") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
