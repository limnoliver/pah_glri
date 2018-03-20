# compare samples, duplicates, and blanks
qa <- make('blanks') %>%
  mutate(unique_id = paste0(State, "-", STAT_ID)) %>%
  filter(!is.na(Parameter)) %>%
  select(unique_id, FIELD_QC_CODE, PARAM_SYNONYM, RESULT, molwt)

dup <- filter(qa, FIELD_QC_CODE != "FB") %>%
  spread(key = FIELD_QC_CODE, value = RESULT) %>%
  mutate(perc_diff = ((SA-DU)/SA)*100)

dup.totals <- group_by(dup, unique_id) %>%
  summarize(DU_total = sum(DU), SA_total = sum(SA)) %>%
  mutate(total_perc_diff = ((SA_total - DU_total)/SA_total)*100) %>%
  mutate(molwt = NA, PARAM_SYNONYM = "Total PAH") %>%
  rename(perc_diff = total_perc_diff, DU = DU_total, SA = SA_total)

dup <- bind_rows(dup, dup.totals) 
dup.means <- group_by(dup, PARAM_SYNONYM) %>%
  summarize(means = mean(SA)) %>%
  arrange(means)

dup$PARAM_SYNONYM <- factor(dup$PARAM_SYNONYM, levels = dup.means$PARAM_SYNONYM)
p <- ggplot(dup, aes(x = reorder(PARAM_SYNONYM, SA), y = perc_diff)) +
  geom_point(aes(color = unique_id), alpha = 0.5) +
  scale_color_discrete(name = "Site ID")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.minor = element_blank()) +
  labs(x = "Compounds, low to high concentration", y = "% diff between sample and duplicate")

ggsave('10_qa_summary/doc/percdiff_by_compound.png', p)
blanks <- filter(qa, FIELD_QC_CODE == "FB")

ggplot(blanks, aes(x = reorder(PARAM_SYNONYM, molwt), y = RESULT)) +
  geom_point(aes(color = unique_id), alpha = 0.5) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.minor = element_blank())

ggplot(qa, aes(x = reorder(PARAM_SYNONYM, molwt), y = RESULT, group = FIELD_QC_CODE)) +
  geom_point(aes(color = unique_id, shape = FIELD_QC_CODE), alpha = 0.8) +
  scale_y_continuous(trans = 'log10')+
  scale_shape_manual(values = c(15, 1, 17)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.minor = element_blank()) +
  labs(x = "Compounds, low to high MW", y = "Concentration (ppb)")


ggplot(blanks, aes())
