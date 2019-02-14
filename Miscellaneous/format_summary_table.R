# this was an attempt to format a table in R
# ended up abandoning and using Excel to format

library(knitr)
library(kableExtra)
library(formattable)
library(htmltools)
library(webshot)

top_sources %>%
  mutate(`PCA top source` = cell_spec(`PCA top source`, 'html', color = 'white',
                                      background = factor(`PCA top source`, unique_top_sources, my_cols)),
         `Profiles top source` = cell_spec(`Profiles top source`, 'html', color = 'white',
                                           background = factor(`Profiles top source`, unique_top_sources, my_cols)),
         `Parent:Alkyl` = formatter("span",
                                    style = ~ style(display = 'inline-block',
                                                    direction = 'ltr',
                                                    'border-radius' = '4px',
                                                    'padding-right' = "2px",
                                                    'background-color' = 'lightgray',
                                                    width = percent(normalize(parent_alkyl))))(`Parent:Alkyl`),
         `HMW:LMW` = formatter("span",
                               style = ~ style(display = 'inline-block',
                                               direction = 'ltr',
                                               'border-radius' = '4px',
                                               'padding-right' = "2px",
                                               'background-color' = 'lightgray',
                                               width = percent(normalize(HMW_LMW))))(`HMW:LMW`)) %>%
  select(sample, `PCA top source`, 'Profiles top source', 'HMW:LMW', 'Parent:Alkyl') %>%
  kable(format = 'html', escape = F, linesep = "") %>%
  kable_styling(bootstrap_options = c('striped', 'condensed'), 
                full_width = FALSE,
                font_size = 12) %>%
  column_spec(1, width = '2cm', bold = TRUE) %>%
  column_spec(2:3, width = '2cm')# %>%
#cat(., file = '12_top_sources_by_site/doc/summary_fig.html')
#as_image(file = '12_top_sources_by_site/doc/summary_fig.png')
#save_kable(file = '12_top_sources_by_site/doc/summary_fig.png')




export_formattable <- function(f, file, width = "60%", height = '90%', 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(ftable, "12_top_sources_by_site/doc/summary_fig.png")


# output list of unique sources
all.top <- unique(c(top_sources$ratio1, top_sources$ratio2, top_sources$ratio3, top_sources$profile, top_sources$pca))
write.csv(all.top, '12_top_sources_by_site/doc/unique_top_sources.csv', row.names = F)
}