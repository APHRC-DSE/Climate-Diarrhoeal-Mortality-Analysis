setwd("D:/APHRC/Hararghe HDSS-Ethiopia")
list.files()

# libraries
library(dplyr)
library(lubridate)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)
library(tidyr)

# load diarr mortality data
Drh <- read.csv("DiarrhoealMortality_Hararghe.csv")

# == select only columns of interest
new_Drh <- Drh %>%
  dplyr::select(DateofDeath, gender, age, age_group) %>%
  dplyr::rename(Date = DateofDeath)

# ==== Fix the "May-14" typo 
new_Drh <- new_Drh %>%
  mutate(age_group = ifelse(age_group == "May-14", "5-14", age_group))
  

# == fix the date
new_Drh <- new_Drh %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"))

# =======================================================================
# LOAD CLIMATE DATA
#======================================================================
clim <- read.csv("Hararghe_Climatedata_2015_2025.csv")

# = filter date to start from 2015 to 2023
clim <- clim %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  filter(date >= as.Date("2015-01-01") &
           date <= as.Date("2023-12-31")) %>%
  rename(Date = date)


# ============================================================
# AGGREGATE MORTALITY TO WEEKLY COUNTS
# ============================================================

weekly_deaths <- new_Drh %>%
  mutate(
    Date = as.Date(Date, format = "%Y-%m-%d"),
    Week = floor_date(Date, unit = "week", week_start = 1)
  ) %>%
  group_by(Week) %>%
  summarise(
    # Total deaths
    deaths_total = n(),
    
    # Deaths by age group
    deaths_under5 = sum(age_group == "Under 5", na.rm = TRUE),
    deaths_5_14   = sum(age_group == "5-14", na.rm = TRUE),
    deaths_15_49  = sum(age_group == "15-49", na.rm = TRUE),
    deaths_50plus = sum(age_group == "50+", na.rm = TRUE),
    
    # Deaths by gender
    deaths_male   = sum(gender == "Male", na.rm = TRUE),
    deaths_female = sum(gender == "Female", na.rm = TRUE),
    
    # Under 5 deaths by gender
    under5_male   = sum(age_group == "Under 5" & gender == "Male", na.rm = TRUE),
    under5_female = sum(age_group == "Under 5" & gender == "Female", na.rm = TRUE),
    
    # Ages 5-14 deaths by gender
    age5_14_male   = sum(age_group == "5-14" & gender == "Male", na.rm = TRUE),
    age5_14_female = sum(age_group == "5-14" & gender == "Female", na.rm = TRUE),
    
    # Ages 15-49 deaths by gender
    age15_49_male   = sum(age_group == "15-49" & gender == "Male", na.rm = TRUE),
    age15_49_female = sum(age_group == "15-49" & gender == "Female", na.rm = TRUE),
    
    # Ages 50+ deaths by gender
    age50plus_male   = sum(age_group == "50+" & gender == "Male", na.rm = TRUE),
    age50plus_female = sum(age_group == "50+" & gender == "Female", na.rm = TRUE),
    
    .groups = "drop"
  )

# Verify weekly mortality counts
glimpse(weekly_deaths)
head(weekly_deaths, 10)


# ============================================================
# AGGREGATE CLIMATE TO WEEKLY
# Rainfall = SUM; temperature = MEAN
# ============================================================

weekly_clim <- clim %>%
  mutate(
    Date = as.Date(Date),
    Week = floor_date(Date, unit = "week", week_start = 1)
  ) %>%
  group_by(Week) %>%
  summarise(
    Rainfall = if (all(is.na(Rainfall))) NA_real_
    else sum(Rainfall, na.rm = TRUE),
    Tmin = if (all(is.na(Tmin))) NA_real_
    else mean(Tmin, na.rm = TRUE),
    Tavg = if (all(is.na(Tavg))) NA_real_
    else mean(Tavg, na.rm = TRUE),
    Tmax = if (all(is.na(Tmax))) NA_real_
    else mean(Tmax, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# CREATE COMPLETE WEEK SEQUENCE — NO GAPS
# ============================================================

full_weeks <- data.frame(
  Week = seq(
    from = min(weekly_clim$Week, na.rm = TRUE),
    to   = max(weekly_clim$Week, na.rm = TRUE),
    by   = "week"
  )
)


# ============================================================
# MERGE WEEKLY MORTALITY WITH COMPLETE WEEK SEQUENCE
# ============================================================

weekly_deaths <- full_weeks %>%
  left_join(weekly_deaths, by = "Week") %>%
  arrange(Week)


# ============================================================
# MERGE MORTALITY AND CLIMATE DATA
# ============================================================

DiarrClimate_weekly <- full_weeks %>%
  left_join(weekly_deaths, by = "Week") %>%
  left_join(weekly_clim, by = "Week") %>%
  arrange(Week)


# ============================================================
# VERIFY FINAL DATASET
# ============================================================

glimpse(DiarrClimate_weekly)
head(DiarrClimate_weekly, 10)

# === Rearrange Variables
DiarrClimate_weekly <- DiarrClimate_weekly %>%
  dplyr::select(Week, Tmin, Tavg, Tmax, Rainfall, everything())

#==Rename Rainfall to Rain
DiarrClimate_weekly <- DiarrClimate_weekly %>%
  dplyr::rename(Rain = Rainfall)

# save weekly Diarr-climate dataset
write.csv(DiarrClimate_weekly, file = "DiarrClimate_Harar_Weekly.csv", row.names = FALSE)


# =========================================================================
# Descriptive statistics of weekly diarrhoeal deaths and climatic variables
#==========================================================================
summary_long <- DiarrClimate_weekly %>%
  dplyr::select(Rain, Tmin, Tavg, Tmax, deaths_total, deaths_under5) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(
    Mean   = round(mean(Value, na.rm = TRUE), 1),
    SD     = round(sd(Value, na.rm = TRUE), 1),
    P25    = round(quantile(Value, 0.25, na.rm = TRUE), 1),
    Median = round(median(Value, na.rm = TRUE), 1),
    P75    = round(quantile(Value, 0.75, na.rm = TRUE), 1),
    Min    = round(min(Value, na.rm = TRUE), 1),
    Max    = round(max(Value, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(summary_long)


# ============================================================
# TIME SERIES PLOTS — WEEKLY DEATHS, RAINFALL AND TEMPERATURE
# ============================================================

# PREPARE DATA
# ============================================================

# Ensure Week is in the correct date format
DiarrClimate_weekly <- DiarrClimate_weekly %>%
  mutate(Week = as.Date(Week))

# ============================================================
# PLOT 1 — WEEKLY ALL-AGE DEATHS
# ============================================================

p1 <- ggplot(
  DiarrClimate_weekly,
  aes(x = Week, y = deaths_total)
) +
  geom_line(
    color = "#E74C3C",
    alpha = 0.7,
    linewidth = 0.5
  ) +
  geom_smooth(
    method = "loess",
    span = 0.1,
    color = "#C0392B",
    linewidth = 1.0,
    se = FALSE
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "A) Weekly Diarrhoeal Deaths, All Ages",
    x = NULL,
    y = "Number of Deaths"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )


# ============================================================
# PLOT 2 — WEEKLY UNDER-FIVE DEATHS
# ============================================================

p_under5 <- ggplot(
  DiarrClimate_weekly,
  aes(x = Week, y = deaths_under5)
) +
  geom_line(
    color = "#8E44AD",
    alpha = 0.7,
    linewidth = 0.5
  ) +
  geom_smooth(
    method = "loess",
    span = 0.1,
    color = "#6C3483",
    linewidth = 1.0,
    se = FALSE
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "B) Weekly Diarrhoeal Deaths, Children Under Five",
    x = NULL,
    y = "Number of Deaths"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )


# ============================================================
# PLOT 3 — WEEKLY TOTAL RAINFALL
# ============================================================

p2 <- ggplot(
  DiarrClimate_weekly,
  aes(x = Week, y = Rain)
) +
  geom_col(
    fill = "#2980B9",
    alpha = 0.7,
    width = 7
  ) +
  geom_smooth(
    method = "loess",
    span = 0.1,
    color = "#1A5276",
    linewidth = 1.0,
    se = FALSE
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "C) Weekly Total Rainfall",
    x = NULL,
    y = "Rainfall (mm)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    plot.margin = margin(5, 10, 0, 10)
  )


# ============================================================
# PLOT 4 — WEEKLY MEAN TEMPERATURE
# Tmin, Tavg and Tmax
# ============================================================

# Reshape temperature data to long format
temp_long <- DiarrClimate_weekly %>%
  dplyr::select(Week, Tmin, Tavg, Tmax) %>%
  pivot_longer(
    cols = c(Tmin, Tavg, Tmax),
    names_to = "Temperature_Type",
    values_to = "Temperature_C"
  ) %>%
  mutate(
    Temperature_Type = factor(
      Temperature_Type,
      levels = c("Tmax", "Tavg", "Tmin")
    )
  )


p3 <- ggplot(
  temp_long,
  aes(
    x = Week,
    y = Temperature_C,
    color = Temperature_Type
  )
) +
  geom_line(
    alpha = 0.6,
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Tmax" = "#E53935",
      "Tavg" = "#1B5E20",
      "Tmin" = "#FF1744"
    ),
    labels = c(
      "Tmax" = "Maximum",
      "Tavg" = "Mean",
      "Tmin" = "Minimum"
    )
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  labs(
    title = "D) Weekly Mean Temperature",
    x = "Year",
    y = "Temperature (°C)",
    color = "Temperature"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.position = "top",
    plot.margin = margin(5, 10, 5, 10)
  )


# ============================================================
# COMBINE ALL FOUR PLOTS — STACKED VERTICALLY
# ============================================================

combined_plot <- p1 / p_under5 / p2 / p3 +
  plot_layout(
    heights = c(1, 1, 1, 1.2)
  ) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 11,
        color = "grey40",
        hjust = 0
      ),
      plot.caption = element_text(
        size = 9,
        color = "grey50",
        hjust = 1
      )
    )
  )


# Display the combined figure
print(combined_plot)


# ============================================================
# SAVE THE COMBINED FIGURE
# ============================================================

# High-resolution TIFF format (600 dpi)
ggsave(
  filename = "Weekly_Deaths_Rainfall_Temperature.tiff",
  plot = combined_plot,
  width = 12,
  height = 16,
  units = "in",
  dpi = 600,
  compression = "lzw"
)















































