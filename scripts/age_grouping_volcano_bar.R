# scripts/age_grouping_volcano_bar.R
# Create volcano plot and barplot of top effect sizes from limma results

library(ggplot2)
library(ggrepel)

LIMMA_RES <- "outputs/limma_results_trend.csv"
OUT_DIR <- "outputs"
TOP_LABEL <- 10
TOPN_BAR <- 20

res <- read.csv(LIMMA_RES, stringsAsFactors = FALSE)
if(!all(c('logFC','P.Value','adj.P.Val','protein','t') %in% colnames(res))) stop('limma results missing required columns')

# compute -log10 p
res$negLogP <- -log10(res$P.Value)
# significance threshold
res$signif <- res$adj.P.Val < 0.05
# fold-change threshold for highlighting (optional)
fc_thr <- 1
res$highlight <- res$signif & abs(res$logFC) >= fc_thr

# volcano plot
volcano_pdf <- file.path(OUT_DIR, 'limma_volcano_pub.pdf')
volcano_png <- file.path(OUT_DIR, 'limma_volcano_pub.png')

p <- ggplot(res, aes(x=logFC, y=negLogP, color=highlight)) +
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=c('grey70','red')) +
  theme_minimal(base_size = 12) +
  xlab('Log2 fold change (trend)') + ylab('-log10(p-value)') +
  ggtitle('Limma trend: Volcano plot') +
  theme(legend.position='none')

# label top hits by adj.P.Val
top_hits <- head(res[order(res$adj.P.Val), ], n = TOP_LABEL)
p <- p + geom_text_repel(data=top_hits, aes(label=protein), size=3, max.overlaps = 20)

# save
pdf(volcano_pdf, width=7, height=5)
print(p)
dev.off()
png(volcano_png, width=7, height=5, units='in', res=300)
print(p)
dev.off()
cat('Wrote volcano:', volcano_pdf, volcano_png, '\n')

# Barplot: top absolute logFC
res$absFC <- abs(res$logFC)
bar_df <- head(res[order(-res$absFC), ], n = TOPN_BAR)
# estimate SE from t-stat (SE = logFC / t) when t != 0
bar_df$SE <- with(bar_df, ifelse(t != 0, logFC / t, NA))
bar_df$protein <- factor(bar_df$protein, levels = rev(bar_df$protein))

bar_pdf <- file.path(OUT_DIR, 'limma_top_effects_bar_pub.pdf')
bar_png <- file.path(OUT_DIR, 'limma_top_effects_bar_pub.png')

bp <- ggplot(bar_df, aes(x=protein, y=logFC, fill=logFC>0)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=logFC - SE, ymax=logFC + SE), width=0.4, na.rm=TRUE) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  xlab('Protein') + ylab('Log2 fold change (trend)') +
  ggtitle(sprintf('Top %d proteins by |logFC|', TOPN_BAR)) +
  scale_fill_manual(values=c('steelblue','tomato')) +
  theme(legend.position='none')

pdf(bar_pdf, width=8, height=6)
print(bp)
dev.off()
png(bar_png, width=8, height=6, units='in', res=300)
print(bp)
dev.off()
cat('Wrote barplots:', bar_pdf, bar_png, '\n')
