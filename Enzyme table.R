library(tidyverse)
library(gt) # Grammar of Table package to produce publishable tables

rm(list=ls())
setwd("~/General Waring/MRes Project")

enzyme.table <- read.csv("Enzyme table.csv", na.strings=c("","NA"))

# convert all "NaN values in DT to NA
enzyme.table <- enzyme.table %>% mutate_all(~ifelse(is.nan(.), NA, .))

# change "inf" values to NA
enzyme.table <- enzyme.table %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

enzyme.table_tbl <- tbl_df(enzyme.table)

enzyme.table_tbl


n = 0
c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
c_container_width = px(800)
c_table_width = px(650)
c_rn = 30
c_save = TRUE
c_format = "PNG"



gt_table <-
  enzyme.table %>% 
  head(c_rn) %>% 
  gt(groupname_col = F) %>% 
  cols_label(
    enzyme = md("Enzyme"), 
    abbreviation = md("Abbreviation"), 
    EC = md("EC"),
    substrate = md("Substrate"))%>% 
  cols_align(
    align = "left",
    columns = c("enzyme","substrate")) %>%
  cols_align(
    align = "center",
    columns = c("abbreviation","EC")) %>% 
  cols_width(
    (enzyme) ~ px(150),
    (abbreviation) ~ px(100),
    (EC) ~ px(100),
    (Function) ~ px(150),
    (substrate) ~ px(150)) %>% 
  tab_options(
    table.font.name = "Optima",
    table.font.color = c_col[1],
    table.border.top.style = "none",
    table.border.bottom.style = "solid",
    table.border.bottom.color = c_col[1],
    table.border.bottom.width = px(3),
    column_labels.border.top.color = "white",
    column_labels.border.top.width = px(3),
    column_labels.border.bottom.color = c_col[1],
    column_labels.border.bottom.width = px(3),
    data_row.padding = px(10))%>%
  tab_style(
    style = list( 
      cell_text(
        size = px(16),
        color = "#2f5375",
        font = "Bloomsbury"
      )
    ),
    locations = list(
      cells_column_labels(everything())
    )
  )  %>% 
  tab_style(
    style = list( 
      cell_text(
        size = px(12),
        color = "#2f5375",
        font = "Bloomsbury"
      ),
      cell_borders(
        sides = c("bottom", "right", "top"),
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
        font = "Bloomsbury", size = px(12), 
        color = "#2f5375")
    ),
    location = list(
      cells_body(columns = c(enzyme))
    )
  ) %>% 
  tab_style(
    style = list(
      cell_text(
        size = px(15),
      ),
      cell_borders(
        sides = c("bottom", "top"),
        color = "white",
        weight = px(1)
      )
    ),
    locations = list(
      cells_body(gt::everything())
    )
  ) 

gt_table

