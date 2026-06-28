################################################################################
# RESPIRATORY VIRUS MORTALITY - FIGURE 1
# UPDATED: H1N1 gap (Nov 2009 - Aug 2010) shown as break in panels B and C
#          Deaths corrected for reporting lag (Jan-Jun 2024)
#          Provisional 2025 deaths shown as dotted line
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
# LOAD DATA
################################################################################

dat <- read.csv("data/ontario_mortality_virology_EXTENDED_1991_2025_CORRECTED.csv",
                stringsAsFactors = FALSE)
dat$date <- as.Date(dat$date)
# Ensure rows exist for the H1N1 gap months so the line breaks (not just shaded)
gap_months <- seq(as.Date("2009-11-01"), as.Date("2010-08-31"), by = "month")
missing <- gap_months[!gap_months %in% dat$date]
if (length(missing) > 0) {
  add <- dat[1:length(missing), ]
  add[,] <- NA
  add$date <- missing
  dat <- rbind(dat, add)
  dat <- dat[order(dat$date), ]
}
# H1N1 exclusion period
h1n1_start <- as.Date("2009-11-01")
h1n1_end   <- as.Date("2010-08-31")

# Set virology to NA for H1N1 gap so ggplot breaks the line
dat <- dat %>%
  mutate(
    monthly_flua_pct_pos = ifelse(date >= h1n1_start & date <= h1n1_end,
                                  NA, monthly_flua_pct_pos),
    monthly_flub_pct_pos = ifelse(date >= h1n1_start & date <= h1n1_end,
                                  NA, monthly_flub_pct_pos),
    monthly_rsv_pct_pos  = ifelse(date >= h1n1_start & date <= h1n1_end,
                                  NA, monthly_rsv_pct_pos)
  )

# Deaths per 100,000
dat <- dat %>%
  mutate(deaths100k = (deaths / (population / 12)) * 100000)

# Key dates
pandemic_start <- as.Date("2020-03-15")

# Shading rectangle for H1N1 gap (used in panels B and C)
h1n1_rect <- data.frame(
  xmin = h1n1_start, xmax = h1n1_end,
  ymin = -Inf,        ymax = Inf
)

################################################################################
# PANEL A: Deaths per 100,000
################################################################################

df_statcan   <- dat %>% filter(lubridate::year(date) < 1994,  !is.na(deaths))
df_odr       <- dat %>% filter(lubridate::year(date) >= 1994, !is.na(deaths),
                                date <= as.Date("2025-02-28"))
df_prov      <- dat %>% filter(date >= as.Date("2025-03-01"),  !is.na(deaths))

# Note: lubridate not always available; use format() as fallback
df_statcan   <- dat %>% filter(as.integer(format(date, "%Y")) < 1994,  !is.na(deaths))
df_odr       <- dat %>% filter(as.integer(format(date, "%Y")) >= 1994, !is.na(deaths),
                                date <= as.Date("2025-02-28"))
df_prov      <- dat %>% filter(date >= as.Date("2025-03-01"),  !is.na(deaths))

p1 <- ggplot() +
  geom_line(data = df_statcan, aes(x = date, y = deaths100k, color = "Statistics Canada"),
            linewidth = 0.7) +
  geom_line(data = df_odr, aes(x = date, y = deaths100k, color = "Ontario Deaths Registry"),
            linewidth = 0.7) +
  geom_line(data = df_prov, aes(x = date, y = deaths100k, color = "Ontario Deaths Registry"),
            linewidth = 0.7, linetype = "dotted", alpha = 0.5) +
  geom_vline(xintercept = pandemic_start, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = c("Statistics Canada"       = "#44AAAA",
                                "Ontario Deaths Registry" = "#E05555")) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(y = "Deaths per 100,000", x = "Date", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

################################################################################
# PANEL B: Influenza A and B (with H1N1 gap)
################################################################################

df_flu <- dat %>%
  select(date, monthly_flua_pct_pos, monthly_flub_pct_pos) %>%
  pivot_longer(cols = c(monthly_flua_pct_pos, monthly_flub_pct_pos),
               names_to = "virus", values_to = "positivity") %>%
  mutate(virus = recode(virus,
                        monthly_flua_pct_pos = "Influenza A",
                        monthly_flub_pct_pos = "Influenza B"))

p2 <- ggplot(df_flu, aes(x = date, y = positivity, color = virus)) +
  # H1N1 gap shading
  geom_rect(data = h1n1_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "gray70", alpha = 0.2) +
  annotate("text", x = as.Date("2010-03-15"), y = Inf,
           label = "H1N1\ngap", vjust = 1.3, size = 2.5, color = "gray50") +
  # Lines — NAs in gap period cause natural break
  geom_line(linewidth = 0.7, na.rm = FALSE) +
  geom_vline(xintercept = pandemic_start, linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = c("Influenza A" = "#E05555",
                                "Influenza B" = "#44AAAA")) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(y = "Percent Positive", x = "Date", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

################################################################################
# PANEL C: RSV (with H1N1 gap)
################################################################################

df_rsv <- dat

p3 <- ggplot(df_rsv, aes(x = date, y = monthly_rsv_pct_pos)) +
  geom_rect(data = h1n1_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "gray70", alpha = 0.2) +
  annotate("text", x = as.Date("2010-03-15"), y = Inf,
           label = "H1N1\ngap", vjust = 1.3, size = 2.5, color = "gray50") +
  geom_line(color = "#E07722", linewidth = 0.7, na.rm = FALSE) +
  geom_vline(xintercept = pandemic_start, linetype = "dashed", linewidth = 0.6) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(y = "RSV Percent Positive", x = "Date") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

################################################################################
# PANEL D: SARS-CoV-2 (dual y-axis)
################################################################################

df_cases <- dat %>%
  filter(date >= as.Date("2020-03-01") & date <= as.Date("2022-08-31") &
           !is.na(monthly_adj_cases) & monthly_adj_cases > 0) %>%
  mutate(adj_cases_10k = monthly_adj_cases / 10000)

df_pos <- dat %>%
  filter(date >= as.Date("2022-09-01") & date <= as.Date("2025-02-28") &
           !is.na(monthly_sars_pct_pos) & monthly_sars_pct_pos > 0)

scale_factor <- max(df_cases$adj_cases_10k, na.rm = TRUE) /
                max(df_pos$monthly_sars_pct_pos, na.rm = TRUE)

p4 <- ggplot() +
  geom_line(data = df_cases,
            aes(x = date, y = adj_cases_10k, color = "Adjusted Cases per 10,000"),
            linewidth = 0.7) +
  geom_line(data = df_pos,
            aes(x = date, y = monthly_sars_pct_pos * scale_factor,
                color = "SARS-2 Positivity (%)"),
            linewidth = 0.7) +
  scale_y_continuous(
    name = "Adjusted Cases per 10,000",
    sec.axis = sec_axis(~ . / scale_factor, name = "SARS-2 Positivity (%)")
  ) +
  scale_color_manual(values = c("Adjusted Cases per 10,000" = "#E05555",
                                "SARS-2 Positivity (%)"     = "#44AAAA")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y",
               guide = guide_axis(angle = 45)) +
  labs(x = "Date", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

################################################################################
# COMBINE AND SAVE
################################################################################

figure1 <- (p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

print(figure1)

ggsave("figures/figure1_mortality_virology.png",  figure1, width = 10, height = 8, dpi = 300)
ggsave("figures/figure1_mortality_virology.tiff", figure1, width = 10, height = 8, dpi = 300)

cat("Figure 1 saved.\n")

sum(dat$date >= as.Date("2009-11-01") & dat$date <= as.Date("2010-08-31"))

