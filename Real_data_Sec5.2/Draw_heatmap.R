source("../EFM algorithm.R")
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)

Y = readRDS("FRED-MD.rds") 

fred_group_num <- c(
  1, 1, 4, 4, 4,              # 1-5
  rep(1, 14),                 # 6-19
  rep(2, 28),                 # 20-47
  rep(3, 10),                 # 48-57
  rep(4, 5),                  # 58-62
  rep(5, 10),                 # 63-72
  rep(8, 3),                  # 73-75
  rep(6, 22),                 # 76-97
  rep(7, 20),                 # 98-117
  rep(2, 3),                  # 118-120
  4,                          # 121
  rep(5, 3),                  # 122-124
  8                           # 125
)

######################################################
Y = Y[,order(fred_group_num, decreasing = FALSE)]
Y_group = sort(fred_group_num, decreasing = FALSE)

set.seed(0)
res50 = efm_est(Y, r=5, tau=0.5, tol=1e-5)

set.seed(0)
res95 = efm_est(Y, r=7, tau=0.95, tol=1e-5)

set.seed(0)
res5 = efm_est(Y, r=6, tau=0.05, tol=1e-5)


# (1) res -> df 
make_df <- function(lmat) {
  mat <- unclass(varimax(lmat)$loadings)
  rownames(mat) <- paste0("R", 1:125)
  colnames(mat) <- paste0("F", 1:ncol(mat))
  
  df <- melt(mat)
  colnames(df) <- c("Row", "Col", "Value")
  
  df$Row    <- factor(df$Row, levels = paste0("R", 1:125))
  df$RowIdx <- as.numeric(df$Row)
  
  a <- quantile(abs(df$Value), 0.8)
  df$WhiteFlag <- abs(df$Value) < a
  df$loading   <- df$Value
  df
}

df_05 <- make_df(res5$lmat)
df_50 <- make_df(res50$lmat)
df_95 <- make_df(res95$lmat)

# set scaling bar and borders
L <- max(abs(c(df_05$loading, df_50$loading, df_95$loading)), na.rm = TRUE)

row_labels <- rep(" ", 125)
row_labels[seq(1, 125, by = 2)] <- paste0("v", seq(1, 125, by = 2))

row_group <- data.frame(
  Row    = row_labels,
  Group  = Y_group,
  RowIdx = 1:125
)

border <- row_group %>%
  group_by(Group) %>%
  summarise(y = max(RowIdx) + 0.5, .groups = "drop")
border_plot <- border[-nrow(border), ]

# plot
make_heatmap_plot <- function(df, title_text) {
  ggplot(df, aes(x = Col, y = RowIdx)) +
    geom_tile(aes(fill = loading), width = 1) +
    scale_y_reverse(
      breaks = row_group$RowIdx,
      labels = row_group$Row
    ) +
    scale_fill_gradientn(
      colours = colorRampPalette(c("green", "white", "red"))(200),
      name = NULL
      #,limits  = c(-L, L)    
    ) +
    geom_tile(
      data = df[df$WhiteFlag, ],
      aes(x = Col, y = RowIdx),
      fill = "white"
    ) +
    geom_hline(
      data = border_plot,
      aes(yintercept = y),
      linewidth = 0.5
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 10),
      axis.title  = element_blank(),
      plot.title = element_text(hjust = 0.5, size=11),
      legend.text = element_text(size = 7),
      legend.position = "right" 
    ) +
    ggtitle(title_text) +
    guides(fill = guide_colorbar(barheight = unit(60, "pt"),
                                 barwidth = unit(10, "pt"))) 
}

p1 <- make_heatmap_plot(df_50, expression("Mean-level loadings (" * alpha == 0.5 * ")"))
p2 <- make_heatmap_plot(df_05, expression("Lower-tail loadings (" * alpha == 0.05 * ")"))
p3 <- make_heatmap_plot(df_95, expression("Upper-tail loadings (" * alpha == 0.95 * ")"))

# three plots
plot_grid(p1, p2, p3, nrow = 1, align = "h")

