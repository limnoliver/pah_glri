
make_map <- function(filename, sample_info, thresh_dat) {
  
  get_flowlines <- function(streamorder, mapRange){
    postURL <- "https://cida.usgs.gov/nwc/geoserver/nhdplus/ows"
    
    filterXML <- paste0('<?xml version="1.0"?>',
                        '<wfs:GetFeature xmlns:wfs="http://www.opengis.net/wfs" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:gml="http://www.opengis.net/gml" service="WFS" version="1.1.0" outputFormat="shape-zip" xsi:schemaLocation="http://www.opengis.net/wfs http://schemas.opengis.net/wfs/1.1.0/wfs.xsd">',
                        '<wfs:Query xmlns:feature="https://gov.usgs.cida/nhdplus" typeName="feature:nhdflowline_network" srsName="EPSG:4326">',
                        '<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                        '<ogc:And>',
                        '<ogc:PropertyIsGreaterThan>',
                        '<ogc:PropertyName>streamorde</ogc:PropertyName>',
                        '<ogc:Literal>',streamorder-1,'</ogc:Literal>',
                        '</ogc:PropertyIsGreaterThan>',
                        '<ogc:BBOX>',
                        '<ogc:PropertyName>the_geom</ogc:PropertyName>',
                        '<gml:Envelope>',
                        '<gml:lowerCorner>',mapRange[3]," ",mapRange[1],'</gml:lowerCorner>',
                        '<gml:upperCorner>',mapRange[4]," ",mapRange[2],'</gml:upperCorner>',
                        '</gml:Envelope>',
                        '</ogc:BBOX>',
                        '</ogc:And>',
                        '</ogc:Filter>',
                        '</wfs:Query>',
                        '</wfs:GetFeature>')
    
    destination = file.path(tempdir(),"nhdflowline_network.zip")
    file <- POST(postURL, body = filterXML, write_disk(destination, overwrite=T))
    
    filePath <- tempdir()
    print("unzipping...")
    unzip(destination, exdir = filePath)
    
    flowLines <- st_read(filePath, layer = 'nhdflowline_network')
    
    return(flowLines)
  }
  
  getBasin <- function(sites, filePath = NA){
    
    postURL <- "https://cida.usgs.gov/nwc/geoserver/NWC/ows"
    
    filterXML <- paste0('<?xml version="1.0"?>',
                        '<wfs:GetFeature xmlns:wfs="http://www.opengis.net/wfs" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:gml="http://www.opengis.net/gml" service="WFS" version="1.1.0" outputFormat="shape-zip" xsi:schemaLocation="http://www.opengis.net/wfs http://schemas.opengis.net/wfs/1.1.0/wfs.xsd">',
                        '<wfs:Query xmlns:feature="https://owi.usgs.gov/NWC" typeName="feature:epa_basins" srsName="EPSG:4326">')
    
    
    if(length(sites) > 1){
      siteText <- ""
      for(i in sites){
        siteText <- paste0(siteText,'<ogc:PropertyIsEqualTo  matchCase="true">',
                           '<ogc:PropertyName>site_no</ogc:PropertyName>',
                           '<ogc:Literal>',i,'</ogc:Literal>',
                           '</ogc:PropertyIsEqualTo>')
      }
      
      filterXML <- paste0(filterXML,'<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                          '<ogc:Or>',siteText,'</ogc:Or>',
                          '</ogc:Filter>')
      
    } else {
      filterXML <- paste0(filterXML,
                          '<ogc:Filter xmlns:ogc="http://www.opengis.net/ogc">',
                          '<ogc:PropertyIsEqualTo matchCase="true">',
                          '<ogc:PropertyName>site_no</ogc:PropertyName>',
                          '<ogc:Literal>',sites,'</ogc:Literal>',
                          '</ogc:PropertyIsEqualTo>',
                          '</ogc:Filter>')
    }
    
    filterXML <- paste0(filterXML,'</wfs:Query>',
                        '</wfs:GetFeature>')
    
    destination = tempfile(pattern = 'basins_shape', fileext='.zip')
    
    file <- POST(postURL, body = filterXML, write_disk(destination, overwrite=T))
    if(is.na(filePath)){
      filePath <- tempdir()
    }
    
    unzip(destination, exdir = filePath)
    basins = st_read(filePath, layer='epa_basins')
    return(basins)
  }
  
  getLakes <- function(){
    shapefile_loc <- "https://www.glahf.org/download/159/"
    
    destination = file.path(tempdir(),"glin_gl_mainlakes.zip")
    file <- GET(shapefile_loc, write_disk(destination, overwrite=T))
    filePath <- tempdir()
    unzip(destination, exdir = filePath)
    
    lakes <- st_read(file.path(filePath, "shoreline_delineations.gdb"), layer = 'shoreline')
    return(lakes)
  }
  
  # set some map features
  mapRange <- c(-93.298145, -74.816895,40.937378, 48.058570 )
  streamorder <- 5
  crsLONGLAT <- 4326
  crs_plot <- st_crs(102003)
  
  # get rivers in map region
  flowlines <- get_flowlines(streamorder, mapRange)
  
  # get great lakes polygons
  lakes <- getLakes()
  
  # set the state and county names of interest
  state_names <- c("minnesota","wisconsin", "michigan","ohio",
                   "illinois","pennsylvania","new york","indiana")
  
  # get STATE data
  GL <- us_states(resolution = "high", states = state_names) %>%
    st_transform(crs = crsLONGLAT)
  
  bb <- st_sfc(
    st_point(mapRange[c(1,3)]),
    st_point(mapRange[c(2,4)]),
    crs = crsLONGLAT) 
  
  bb_proj <- st_transform(bb, crs = crs_plot)
  b <- st_bbox(bb_proj)
  
  # create categories for color plotting
  thresh_dat <- thresh_dat[[1]]
  
  # remove sites that don't have all thresholds
  # (this is the case where esbtu was NA when there was no
  # TOC data)
  
  
  thresh_dat$cat <- NA
  for (i in 1:nrow(thresh_dat)) {
    
    # skip rows with NA vals in thresholds
    if (complete.cases(thresh_dat[i, c('tec_ratio', 'pec_ratio', 'sum_esbtu')]) == FALSE) {
      next
    }
    
  if (thresh_dat$tec_ratio[i] < 1 & thresh_dat$sum_esbtu[i] < 1) {
    thresh_dat$cat[i] <- '< all thresholds'
  } else if (thresh_dat$tec_ratio[i] > 1 & thresh_dat$pec_ratio[i] < 1 & thresh_dat$sum_esbtu[i] < 1){
    thresh_dat$cat[i] <- '> TEC'
  } else if (thresh_dat$pec_ratio[i] > 1 & thresh_dat$sum_esbtu[i] < 1){
    thresh_dat$cat[i] <- '> PEC'
  } else if (thresh_dat$pec_ratio[i] > 1 & thresh_dat$sum_esbtu[i] > 1){
    thresh_dat$cat[i] <- '> all thresholds'
  } else if (thresh_dat$sum_esbtu[i] > 1 & thresh_dat$pec_ratio[i] <1 & thresh_dat$tec_ratio[i] <1) {
    thresh_dat$cat[i] <- '> ESBTU'
  } else if (thresh_dat$sum_esbtu[i] > 1 & thresh_dat$pec_ratio[i] < 1 & thresh_dat$tec_ratio[i] >1) {
    thresh_dat$cat[i] <- '> ESBTU & TEC'
  }
  }
  
  thresh_dat <- filter(thresh_dat, !is.na(cat))
  
  # create sites data that will be used for plotting
  sites_df <- left_join(thresh_dat, distinct(sample_info[,c('unique_id','Lat', 'Lon')]), by = "unique_id")
  
  sites_df <- st_as_sf(sites_df,
                       coords = c("Lon","Lat"),
                       crs = crsLONGLAT)
  sites_proj <- st_transform(sites_df, crs = crs_plot)
  
  my_cols <- brewer.pal(n = 4, 'YlOrRd')
  
  my_cols <- c('#808080', my_cols[c(4,2,3,2,1)])
  # make map
  p <- ggplot() + 
    geom_sf(data = lakes, fill = "lightblue", color = "lightblue") +
    geom_sf(data = flowlines, color = "lightblue") +
    geom_sf(data = GL, color = "gray50", fill=NA) +
    geom_sf(data = sites_proj, alpha = 0.4, shape = 21, 
            aes(size = sum_EPA16, fill = cat), show.legend = FALSE) +
    scale_size(range = c(3,14), guide = FALSE, breaks = c(500, 15000, 150000)) +
    scale_fill_manual(values = my_cols) + 
    coord_sf(crs = crs_plot,
             xlim = c(b["xmin"],b["xmax"]), 
             ylim = c(b["ymin"],b["ymax"])) +
    theme_minimal() +
    theme(panel.grid.major = element_line(colour = 'transparent'), #Bug that's apparently getting fixed
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(), 
          #legend.position = c(.75, .75),
          legend.direction = 'horizontal') +
    geom_point(alpha = 0, shape = 16,
               aes(x = rep(214731.983109861, nrow(thresh_dat)), y = rep(838589.951769598, nrow(thresh_dat)), 
                   size = sites_proj$sum_EPA16, fill = sites_proj$cat)) + # this is a hack to fix interaction between sf and ggplot2 legend bug
    guides(size = guide_legend(label.position = 'bottom', label.hjust = 0.5,  
                               title = 'Sum PAH16 (ppb)', title.position = 'top', 
                               override.aes = list(alpha = 1, stroke = 2), order = 2), 
           fill = guide_legend(title.position = 'top', title = "Sample Relative to Thresholds", order = 1, ncol = 2,
                               override.aes = list(alpha = 0.7, size = 4, color = my_cols))) + 
    theme(legend.box = "vertical")
  
  ggsave(filename, p, height = 6, width = 10)
}


