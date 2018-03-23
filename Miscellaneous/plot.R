head(samples)
library(ggplot2)
library(dplyr)
library(pah)



ggplot(df, aes(x = reorder(unique_id, RESULT), RESULT)) +
  geom_bar(stat="identity", position="identity", colour="black", aes(fill = State)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~PARAM_SYNONYM, ncol = ifelse(length(unique(df$PARAM_SYNONYM)) > 3, 2, 1),
             scales = 'free')
           

p <- ggplot(df, aes(x = reorder(unique_id, RESULT), y = RESULT)) +
  geom_bar(stat="identity", position="identity", colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p$mapping$fill <- df[["State"]]
p
p + aes(fill = State)

group_columns <- c("State", "Watershed")

df <- df %>%
  ungroup() %>%
  arrange(Lake, Watershed, RESULT) %>%
  mutate(unique_id = factor(unique_id, levels = unique(unique_id)))

state.order <- df %>% 
                group_by(Lake, Watershed) %>% 
                summarize(mean = mean(RESULT)) %>%
                arrange(mean)
df$Watershed <- factor(df$Watershed, levels = state.order$Watershed)

tick_labeller <- function(x) {
  lab <- gsub("MN-|OH-|MI-|WI-|NY-|IN-", "", x)
}

df$RESULT <- df$RESULT*1000 # modify from ng/g to mg/kg
  
df <- df %>%
  mutate(threshold = ifelse(PARAM_SYNONYM == "Priority Pollutant PAH", log10(22800), NA))

dummy.dat <- df %>%
  select(State, PARAM_SYNONYM, threshold) %>%
  distinct()

p <- ggplot(df, aes(x = reorder(unique_id, RESULT), y = log10(RESULT), fill = Lake)) +
  geom_bar(stat="identity", width = 0.9, position = position_dodge(width = 2), colour="black", 
           aes(fill = Lake)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #facet_wrap(~PARAM_SYNONYM, ncol = ifelse(length(unique(df$PARAM_SYNONYM)) > 3, 2, 1),
  #           scales = 'free') +
  facet_grid(.~Watershed, space = "free_x", scales = 'free_x', switch = "x") +
  # probably have to add below line in glri workflow
  geom_hline(data = dummy.dat, aes(yintercept = threshold)) +
  
  theme(panel.spacing=unit(5,"pt"), strip.background = element_blank(),
        strip.text = element_text(margin = margin(4,10,4,10, unit="pt")), 
        axis.text.x = element_text(margin = margin(t = -120, b = 100)), axis.ticks.x = element_blank(), 
        strip.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Concentration (ppb)") +
  # will need to do below line in glri workflow
  scale_x_discrete(labels=tick_labeller) +
  scale_y_continuous(limits = c(0, 5),breaks = c(1,2,3,4,5), 
                     labels = c(10, 100,1000, 10000,100000))

ggsave("test.pdf", p, height = 6, width = 14)
