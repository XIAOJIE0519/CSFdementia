library(ggplot2)
library(dplyr)
library(tidyr)

df <- read.csv('./compare/incremental_auc_data.csv', stringsAsFactors = FALSE)

bio_order <- c('pTau', 'tTau', 'pTau181', 'pTau217', 'AB42', 'AB40')
study_colors <- c(
  study_1  = '#1f77b4',
  study_4  = '#ff7f0e',
  study_9  = '#d62728',
  study_11 = '#9467bd'
)

df$Biomarker <- factor(df$Biomarker, levels = bio_order)
df$bio_base  <- as.numeric(df$Biomarker)

df_long <- df |>
  pivot_longer(cols = c(AUC_Base, AUC_Plus_EODstage),
               names_to = 'Type', values_to = 'AUC') |>
  mutate(
    Type = ifelse(Type == 'AUC_Base', 'Base', '+EODstage'),
    Type = factor(Type, levels = c('Base', '+EODstage'))
  )

# y-positions: each biomarker block, studies offset within
df_long <- df_long |>
  group_by(Biomarker, Study) |>
  mutate(bio_base = as.numeric(Biomarker)) |>
  ungroup()

# Background rectangles
bio_bg <- data.frame(
  ymin = seq_along(bio_order) - 0.45,
  ymax = seq_along(bio_order) + 0.45,
  fill = rep(c('#f5f5f5', '#ffffff'), length.out = length(bio_order))
)

p <- ggplot() +
  geom_rect(
    data = bio_bg,
    aes(xmin = 0.38, xmax = 1.03, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.4, inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_segment(
    data = df,
    aes(x = AUC_Base, xend = AUC_Plus_EODstage,
        y = bio_base, yend = bio_base, color = Study),
    linewidth = 1.4, alpha = 0.75
  ) +
  geom_point(
    data = df_long,
    aes(x = AUC, y = bio_base, color = Study, shape = Type),
    size = 4.5, alpha = 0.95
  ) +
  scale_color_manual(values = study_colors) +
  scale_shape_manual(values = c('Base' = 16, '+EODstage' = 17)) +
  scale_x_continuous(limits = c(0.38, 1.03),
                     breaks = seq(0.4, 1.0, 0.1),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(bio_order),
                     labels = bio_order,
                     expand = c(0.04, 0.04)) +
  labs(
    x     = 'AUC',
    y     = '',
    title = 'Incremental AUC: Biomarker Alone vs Biomarker + EODstage',
    subtitle = 'Circle = Biomarker alone  |  Triangle = + EODstage (LogisticRegression)',
    color = 'Study',
    shape = 'Condition'
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(face = 'bold', size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = 'grey40'),
    axis.title.x    = element_text(face = 'bold', size = 11),
    axis.text.y     = element_text(face = 'bold', size = 11),
    axis.text.x     = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = 'right',
    legend.title       = element_text(face = 'bold'),
    plot.margin        = margin(10, 15, 10, 10)
  )

ggsave('./compare/incremental_lollipop.png', p,
       width = 11, height = 7, dpi = 300)
cat('Saved: incremental_lollipop.png\n')

ggsave('./compare/incremental_lollipop.pdf', p,
       width = 11, height = 7)
cat('Saved: incremental_lollipop.pdf\n')
