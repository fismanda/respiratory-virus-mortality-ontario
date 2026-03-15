# ============================================================
# Figure 3 – Forest plot: PAF for respiratory viruses
# Fisman et al. Population Attributable Mortality, Ontario
# ============================================================
# Requires: ggplot2, dplyr
# install.packages(c("ggplot2", "dplyr"))

library(ggplot2)
library(dplyr)

# ── 1. DerSimonian-Laird random-effects meta-analysis ───────
dl_meta <- function(par, lcl, ucl) {
  se   <- (ucl - lcl) / 3.919928
  vi   <- se^2
  wi   <- 1 / vi
  fe   <- sum(wi * par) / sum(wi)
  Q    <- sum(wi * (par - fe)^2)
  k    <- length(par)
  cc   <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- max(0, (Q - (k - 1)) / cc)
  wr   <- 1 / (vi + tau2)
  pool <- sum(wr * par) / sum(wr)
  se_p <- sqrt(1 / sum(wr))
  I2   <- if (Q > 0) max(0, (Q - (k - 1)) / Q * 100) else 0
  list(pool = pool,
       lcl  = pool - 1.959964 * se_p,
       ucl  = pool + 1.959964 * se_p,
       I2   = I2,
       wi   = wr)
}

# ── 2. Colour palette ────────────────────────────────────────
pal <- c(
  pre_nofft  = "#1a56db",
  pre_fft    = "#76a9fa",
  pan_nofft  = "#b45309",
  pan_fft    = "#d97706",
  post_nofft = "#9f1239",
  post_fft   = "#e11d48",
  pooled     = "#1e1e1e"
)

# ── 3. Data ──────────────────────────────────────────────────
raw <- tribble(
  ~section,                       ~label,                    ~par,      ~lcl,      ~ucl,      ~colour,
  "All Respiratory Viruses",      "Pre-pandemic, no FFT",    0.048727,  0.043367,  0.054057,  "pre_nofft",
  "All Respiratory Viruses",      "Pre-pandemic, with FFT",  0.015305,  0.004546,  0.025948,  "pre_fft",
  "All Respiratory Viruses",      "Pandemic, no FFT",        0.108060,  0.084804,  0.130725,  "pan_nofft",
  "All Respiratory Viruses",      "Pandemic, with FFT",      0.076115,  0.052288,  0.099344,  "pan_fft",
  "Influenza A",                  "Pre-pandemic, no FFT",    0.028892,  0.024864,  0.032903,  "pre_nofft",
  "Influenza A",                  "Pre-pandemic, with FFT",  0.018211,  0.013891,  0.022512,  "pre_fft",
  "Influenza A",                  "Pandemic, no FFT",        0.011574,  0.004111,  0.018981,  "pan_nofft",
  "Influenza A",                  "Pandemic, with FFT",      0.006349, -0.000110,  0.012767,  "pan_fft",
  "Influenza B",                  "Pre-pandemic, no FFT",    0.002131, -0.000692,  0.004946,  "pre_nofft",
  "Influenza B",                  "Pre-pandemic, with FFT",  0.001533, -0.001355,  0.004411,  "pre_fft",
  "Influenza B",                  "Pandemic, no FFT",        0.003514, -0.002578,  0.009568,  "pan_nofft",
  "Influenza B",                  "Pandemic, with FFT",      0.000087, -0.005621,  0.005763,  "pan_fft",
  "Respiratory Syncytial Virus",  "Pre-pandemic, no FFT",    0.018650,  0.013308,  0.023962,  "pre_nofft",
  "Respiratory Syncytial Virus",  "Pre-pandemic, with FFT", -0.004529, -0.012104,  0.002990,  "pre_fft",
  "Respiratory Syncytial Virus",  "Pandemic, no FFT",        0.020157,  0.010622,  0.029599,  "pan_nofft",
  "Respiratory Syncytial Virus",  "Pandemic, with FFT",      0.009467, -0.000135,  0.018977,  "pan_fft",
  "SARS-CoV-2",                   "PHEIC, no FFT",           0.093000,  0.069000,  0.116000,  "pan_nofft",
  "SARS-CoV-2",                   "PHEIC, with FFT",         0.073000,  0.050000,  0.096000,  "pan_fft",
  "SARS-CoV-2",                   "Post-PHEIC, no FFT",      0.083000,  0.024000,  0.139000,  "post_nofft",
  "SARS-CoV-2",                   "Post-PHEIC, with FFT",    0.098000,  0.011000,  0.177000,  "post_fft"
)

# ── 4. Pooled estimates + RE weights ─────────────────────────
section_order <- c("All Respiratory Viruses", "Influenza A",
                   "Influenza B", "Respiratory Syncytial Virus", "SARS-CoV-2")

raw$weight  <- NA_real_
raw$is_pool <- FALSE
raw$is_title <- FALSE
pooled_rows <- list()

for (sec in section_order) {
  idx <- which(raw$section == sec)
  m   <- dl_meta(raw$par[idx], raw$lcl[idx], raw$ucl[idx])
  raw$weight[idx] <- m$wi
  exp_flag <- if (sec == "SARS-CoV-2") "\u2020" else ""
  pooled_rows[[sec]] <- tibble(
    section  = sec,
    label    = sprintf("Pooled (I\u00B2=%.1f%%)%s", m$I2, exp_flag),
    par      = m$pool, lcl = m$lcl, ucl = m$ucl,
    colour   = "pooled", weight = NA_real_,
    is_pool  = TRUE, is_title = FALSE
  )
}
pooled_df <- bind_rows(pooled_rows)

# ── 5. Build y layout (1 unit per row, 0.6 spacer between sections) ──
SPACER <- 0.6
layout_rows <- list()
y <- 0

for (sec in section_order) {
  # title
  y <- y + 1
  layout_rows[[length(layout_rows)+1]] <- tibble(
    y=y, section=sec, label=sec, par=NA, lcl=NA, ucl=NA,
    colour="pooled", weight=NA, is_pool=FALSE, is_title=TRUE
  )
  # data rows
  for (i in seq_len(sum(raw$section == sec))) {
    y <- y + 1
    layout_rows[[length(layout_rows)+1]] <-
      raw[raw$section == sec, ][i, ] %>% mutate(y=y, is_title=FALSE)
  }
  # pooled
  y <- y + 1
  layout_rows[[length(layout_rows)+1]] <-
    pooled_df[pooled_df$section == sec, ] %>% mutate(y=y, is_title=FALSE)
  if (sec != tail(section_order, 1)) y <- y + SPACER
}

layout <- bind_rows(layout_rows)
layout$colour_hex <- pal[layout$colour]
layout$colour_hex[is.na(layout$colour_hex)] <- "#1e1e1e"

# Fixed marker size for all data rows (CI width communicates precision)
layout$pt_size <- 2.5

# PAF label strings
layout <- layout %>% mutate(
  ci_label = if_else(
    !is_title & !is.na(par),
    sprintf("%.1f%% (%.1f%%, %.1f%%)", par*100, lcl*100, ucl*100),
    ""
  )
)

y_max <- max(layout$y) + 0.8
X_MIN <- -0.025; X_MAX <- 0.205

# ── 6. Diamond polygons for pooled rows ──────────────────────
diamond_df <- layout %>%
  filter(is_pool) %>%
  rowwise() %>%
  reframe(
    xd = c(max(lcl, X_MIN), par, min(ucl, X_MAX), par),
    yd = c(y,               y - 0.20, y,           y + 0.20),
    grp = section
  )

# Separator lines above pooled rows
sep_df <- layout %>% filter(is_pool) %>% mutate(y_sep = y - 0.45)

# ── 7. Plot ──────────────────────────────────────────────────
p <- ggplot() +

  # vertical null line
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#999999", linewidth = 0.5) +

  # light x grid
  geom_vline(xintercept = seq(-0.02, 0.20, 0.02),
             colour = "#eeeeee", linewidth = 0.35) +

  # separator above pooled rows
  geom_segment(data = sep_df,
               aes(x = X_MIN, xend = X_MAX, y = y_sep, yend = y_sep),
               colour = "#cccccc", linewidth = 0.4) +

  # CI lines
  geom_segment(
    data = layout %>% filter(!is_pool, !is_title, !is.na(par)),
    aes(x     = pmax(lcl, X_MIN), xend = pmin(ucl, X_MAX),
        y     = y,                yend = y,
        colour = I(colour_hex)),
    linewidth = 0.9, lineend = "round") +

  # Square markers — size in mm, NOT mapped through scale
  geom_point(
    data = layout %>% filter(!is_pool, !is_title, !is.na(par)),
    aes(x = par, y = y, colour = I(colour_hex), size = pt_size),
    shape = 15, show.legend = FALSE) +

  scale_size_identity() +   # treat size values as literal mm

  # Pooled diamonds
  geom_polygon(data = diamond_df,
               aes(x = xd, y = yd, group = grp),
               fill = "#1e1e1e", colour = "#1e1e1e",
               linewidth = 0.4, alpha = 0.88) +

  # Row labels (left of plot, right-aligned)
  geom_text(data = layout %>% filter(!is_title),
            aes(x = X_MIN - 0.003, y = y, label = label,
                fontface = if_else(is_pool, "italic", "plain")),
            hjust = 1, size = 2.65, colour = "#2a2a2a") +

  # Section title labels (bold)
  geom_text(data = layout %>% filter(is_title),
            aes(x = X_MIN - 0.003, y = y, label = label),
            hjust = 1, size = 3.3, fontface = "bold",
            colour = "#111111") +

  # PAF text (right of plot)
  geom_text(data = layout %>% filter(!is_title, !is.na(par)),
            aes(x = X_MAX + 0.003, y = y, label = ci_label,
                fontface = if_else(is_pool, "bold", "plain")),
            hjust = 0, size = 2.55, colour = "#222222",
            family = "mono") +

  # Axes
  scale_x_continuous(
    limits = c(X_MIN - 0.045, X_MAX + 0.055),
    breaks = seq(-0.02, 0.20, 0.02),
    labels = function(x) paste0(round(x * 100), "%"),
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    limits = c(y_max, 0.2),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +

  labs(
    x       = "Population Attributable Fraction (%)",
    caption = paste0(
      "\u2020 SARS-CoV-2 sub-period estimates are exploratory ",
      "(PHEIC n=38, post-PHEIC n=22 months). ",
      "Primary inference rests on the combined pandemic model."
    )
  ) +

  theme_classic(base_size = 9) +
  theme(
    # backgrounds explicitly white
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    # y axis — hide everything
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.title.y       = element_blank(),
    # x axis
    axis.line.x        = element_line(linewidth = 0.5, colour = "#555555"),
    axis.ticks.x       = element_line(linewidth = 0.4, colour = "#555555"),
    axis.text.x        = element_text(size = 7.5, colour = "#333333"),
    axis.title.x       = element_text(size = 8.5, colour = "#111111",
                                      margin = margin(t = 6)),
    # no grid (we drew our own vlines)
    panel.grid         = element_blank(),
    # generous left margin for labels
    plot.margin        = margin(t = 12, r = 8, b = 10, l = 155),
    legend.position    = "none",
    plot.caption       = element_text(size = 6.8, colour = "#555555",
                                      hjust = 0, face = "italic",
                                      margin = margin(t = 8))
  )

# ── 8. Save ──────────────────────────────────────────────────
ggsave("figure3_forest_plot.png", plot = p,
       width = 12, height = 10.5, units = "in",
       dpi = 300, bg = "white")

message("Done — figure3_forest_plot.png saved.")
