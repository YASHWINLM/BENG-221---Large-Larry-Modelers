## ===============================
## Packages
## ===============================
## If not installed, run:
## install.packages("ggplot2")
## install.packages("ggnewscale")

library(ggplot2)
library(ggnewscale)

## ===============================
## 1. 6 taxa × 7 KO (0/1 matrix)
## ===============================

ko_df <- data.frame(
  feature = c("Prevotella",
              "Unclassified S24-7",
              "Unclassified Enterobacteriaceae",
              "Clostridium",
              "Unclassified Lachnospiraceae",
              "Unclassified Christensenellaceae"),
  K00625 = c(1, 0, 0, 1, 1, 1),  # Acetate: pta
  K00925 = c(1, 1, 1, 1, 1, 1),  # Acetate: ackA
  K13922 = c(0, 0, 1, 0, 1, 0),  # Propionate: pduP
  K00929 = c(0, 0, 0, 1, 0, 0),  # Butyrate: buk
  K01034 = c(0, 0, 1, 0, 0, 0),  # Butyrate: CoA-transferase (but/atoD-like)
  K00248 = c(0, 0, 0, 1, 1, 1),  # Butyrate: Bcd
  K00016 = c(0, 0, 0, 1, 1, 1)   # Butyrate: Hbd-like dehydrogenase
)

## ===============================
## 2. Total SCFA coefficients (–1 to 1)
## ===============================

coef_vec <- c(
  "Prevotella"                       =  0.947945665,
  "Unclassified S24-7"               =  0.436813924,
  "Unclassified Enterobacteriaceae"  =  0.837629475,
  "Clostridium"                      =  0.119939472,
  "Unclassified Lachnospiraceae"     =  0.104173274,
  "Unclassified Christensenellaceae" = -0.666180558
)

## y-axis from top to bottom = coef from low to high
## In ggplot, the first factor level is at the bottom,
## so we order levels from high to low
feature_order <- names(coef_vec)[order(coef_vec, decreasing = TRUE)]

## ===============================
## 3. Wide → long (x = KO, y = feature)
## ===============================

ko_mat <- as.matrix(ko_df[, -1])   # remove "feature" column
rownames(ko_mat) <- ko_df$feature

df_long <- as.data.frame(as.table(ko_mat))
colnames(df_long) <- c("feature", "KO", "present")
df_long$present <- as.numeric(as.character(df_long$present))

## Assign each KO to an SCFA pathway (for tile colors)
df_long$SCFA <- ifelse(df_long$KO %in% c("K00625","K00925"), "Acetate",
                ifelse(df_long$KO %in% c("K13922"), "Propionate",
                ifelse(df_long$KO %in% c("K00929","K01034","K00248","K00016"),
                       "Butyrate", NA)))

## KO order: last column will be Total_coef
df_long$KO <- factor(
  df_long$KO,
  levels = c("K00625","K00925",           # Acetate
             "K13922",                    # Propionate
             "K00929","K01034",
             "K00248","K00016",           # Butyrate
             "Total_coef")                # coef column (added below)
)

## feature order: by SCFA coefficient
df_long$feature <- factor(
  df_long$feature,
  levels = feature_order
)

## ===============================
## 4. Data frame for coef column (KO = "Total_coef")
## ===============================

coef_df <- data.frame(
  feature = factor(feature_order, levels = feature_order),
  KO      = factor("Total_coef", levels = levels(df_long$KO)),
  coef    = coef_vec[feature_order]
)

## ===============================
## 5. Plot: no grid, no outer border
## ===============================

p <- ggplot() +
  ## (1) SCFA gene presence tiles
  geom_tile(
    data = subset(df_long, present == 1),
    aes(x = KO, y = feature, fill = SCFA),
    color = NA,            # no tile borders
    linewidth = 0.3
  ) +
  scale_fill_manual(
    name   = "SCFA genes",
    values = c(
      "Acetate"    = "#9467bd",  # purple
      "Propionate" = "#2ca02c",  # green
      "Butyrate"   = "#ff7f0e"   # orange
    ),
    na.translate = FALSE
  ) +
  ## (2) New fill scale for SCFA coefficients
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = coef_df,
    aes(x = KO, y = feature, fill = coef),
    color = NA,
    linewidth = 0.3
  ) +
  scale_fill_gradient2(
    name     = "Total SCFA coef",
    low      = "#2166ac",  # blue
    mid      = "white",
    high     = "#b2182b",  # red
    midpoint = 0,
    limits   = c(-1, 1)
  ) +
  ## Axes and theme
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "KO", y = "Taxon") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),                       # no background grid
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA), # no outer border
    panel.border     = element_blank(),
    plot.margin      = margin(5, 5, 5, 5),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y      = element_text(color = "black"),
    legend.box       = "vertical",
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10)
  )

print(p)

