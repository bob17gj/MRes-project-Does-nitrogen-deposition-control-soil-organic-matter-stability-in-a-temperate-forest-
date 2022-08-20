library(tidyverse)
library(gt) # Grammar of Table package to produce publishable tables

rm(list=ls())
setwd("~/General Waring/MRes Project")

site.properties <- read.csv("Site properties.csv", na.strings=c("","NA"))


# convert all "NaN values in DT to NA
site.properties <- site.properties %>% mutate_all(~ifelse(is.nan(.), NA, .))

# change "inf" values to NA
site.properties <- site.properties %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

# make a "tibble" copy of site.properties

site.properties_tbl <- tbl_df(site.properties)

site.properties


n = 0
c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
c_container_width = px(800)
c_table_width = px(650)
c_rn = 30
c_save = TRUE
c_format = "PNG"


# find mean values

mean.site.properties <- site.properties %>% 
  group_by(site) %>% 
  summarise(mean.GSM = mean(GSM),
            mean.pH = mean(pH.H2O),
            mean.TN = mean(TN),
            mean.TC = mean(TC),
            mean.TOC = mean(TOC),
            mean.N.deposition = mean(NOy.NHx.combined.28.year),
            mean.MAP = mean(MAP),
            mean.MAT = mean(MAT),
            parent.material,
            soil.type,
            forest.type,
            longitude,
            latitude) 

mean.site.properties <- distinct(mean.site.properties) # remove duplicate rows


gt_table <-
  mean.site.properties %>% 
  head(c_rn) %>% 
  gt(groupname_col = F) %>% 
  cols_label(
             site = md("Site"), 
             longitude = md("Longitude"), 
             latitude = md("Latitude"),
             soil.type = md("Soil type"), 
             forest.type = md("Forest type"), 
             parent.material = md("Parent material"), 
             mean.N.deposition = html("N depositon (kg ha<sup>-1</sup>yr<sup>-1</sup>)"),
             mean.TC = md("TC (%)"), 
             mean.TOC = md("TOC (%)"), 
             mean.TN = md("TN (%)"), 
             mean.pH = md("pH (%)"), 
             mean.MAT = md("MAT (\u00B0C)"), 
             mean.MAP = md("MAP (mm)"),
             mean.GSM = md("GSM (%)")) %>% 
  fmt_number(
    columns = c("mean.N.deposition", "mean.TC", "mean.TOC","mean.TN","mean.pH",
                "mean.MAT", "mean.MAP", "mean.GSM", 
                ),
    suffixing = F,
    decimals = 2) %>%
  cols_align(
    align = "center",
    columns = c("parent.material","mean.N.deposition",
                   "mean.MAP","mean.MAT")) %>% 
  cols_align(
    align = "left",
    columns = c("site","longitude","latitude",
                   "mean.GSM","mean.pH",
                   "mean.TN","mean.TC","mean.TOC")) %>% 
  cols_width(
    vars(site) ~ px(60),
    vars(longitude) ~ px(100),
    vars(latitude) ~ px(100),
    vars(soil.type) ~ px(100),
    vars(forest.type) ~ px(120),
    vars(parent.material) ~ px(170),
    vars(mean.N.deposition) ~ px(240),
    vars(mean.TC) ~ px(100),
    vars(mean.TOC) ~ px(120),
    vars(mean.TN) ~ px(100),
    vars(mean.pH) ~ px(100),
    vars(mean.MAT) ~ px(120),
    vars(mean.MAP) ~ px(120),
    vars(mean.GSM) ~ px(100),
  )%>%
  tab_footnote(
    footnote = "NA denotes that data is not available for the site characteristic.",
    locations = 
      cells_column_labels(columns = c("mean.TN","mean.TC","mean.TOC",
                                       "parent.material"))
  ) %>%
  tab_options(
    table.font.name = "Optima",
    table.font.color = c_col[1],
    table.border.top.style = "none",
    table.border.bottom.style = "solid",
    table.border.bottom.color = c_col[2],
    table.border.bottom.width = px(2),
    column_labels.border.top.color = "white",
    column_labels.border.top.width = px(3),
    column_labels.border.bottom.color = c_col[2],
    column_labels.border.bottom.width = px(3),
    data_row.padding = px(10),
    footnotes.font.size = px(15),
    table.font.size = px(16)
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(20)
      ),
      cell_borders(
        sides = c("bottom", "top"),
        color = c_col[1],
        weight = px(1)
      )
    ),
    locations = list(
      cells_body(gt::everything())
    )
  ) %>%
  tab_style(
    style = list( 
      cell_text(
        size = px(22),
        color = "#2f5375",
        font = "Bloomsbury"
      )
    ),
    locations = list(
      cells_column_labels(everything())
    )
  ) %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(16),
        color = "#2f5375",
        font = "Bloomsbury"
      ),
      cell_borders(
        sides = c("bottom"),
        style = "solid",
        color = c_col[1],
        weight = px(2)
      )
    ),
    locations = list(
      cells_row_groups(gt::everything())
    )
  ) %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(40),
        color = "#2f5375",
        font = "Bloomsbury"
      ),
      cell_borders(
        sides = c("bottom", "right"),
        style = "solid",
        color = "white",
        weight = px(1)
      )
    ),
    locations = list(
      cells_stub(gt::everything()),
      cells_stubhead()
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        font = "Bloomsbury", size = px(16), 
        color = "#2f5375")
    ),
    location = list(
      cells_body(columns = c(site))
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(14),
        weight = "normal",
        align = "left",
        font = "Bloomsbury"
      )
    ),
    locations = list(
      cells_title(groups = "title")
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(16),
        align = "left"
      )
    ),
    locations = list(
      cells_title(groups = "subtitle")
    )
  ) %>% 
  tab_footnote(
    footnote = 
    "TC = total C, 
    TOC = total organic C, 
    TN = total N, 
    GSM = gravimetric soil moisture, 
    MAP = mean annual precipitation, 
    MAT = mean annual temperature,
    UD = unconsolidated deposits, 
    UGD = unconsolidated glacial deposits, 
    CCS = consolidated clastic sedimentary rocks, 
    MR = metamorphic rocks.")%>% 
  tab_style(
    style = list(
      cell_text(
        size = px(17),
      ),
      cell_borders(
        sides = c("bottom", "top"),
        color = c_col[1],
        weight = px(1)
      )
    ),
    locations = list(
      cells_body(gt::everything())
    )
  )  
gt_table


# reorder headers
gt_table <- 
  gt_table %>% 
  cols_move(columns = c(latitude, longitude, forest.type, soil.type, 
                        parent.material, mean.N.deposition,
                        mean.TC, mean.TOC, mean.TN, 
                        mean.pH, mean.GSM,
                        mean.MAP,mean.MAT),
            after = site)
gt_table 

#cols_move(longitude,latitude,soil.type,forest.type,parent.material,
 #                                mean.N.deposition,mean.TC,mean.TOC,mean.TN,mean.pH,
  #                               mean.MAT,mean.MAP,mean.GSM,mean.BD,
   #       after = site)


gtsave(gt_table, "Site_properties_6.pdf")

