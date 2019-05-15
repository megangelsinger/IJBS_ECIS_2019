load("gel_data.Rdata")
load("bsa_data.Rdata")
source("ECISfunctions.R")
source("knn_check.R")
source("eor_vs_indiv_tis_check.R")

################################# Analysis ###################################
library(tree)  # for 'tree'
library(MASS) # for 'lda' and 'qda'
library(parallel) # for 'mclapply'

# 9.5 hours to run
results_gel <- mainAnalysis(gel.data)

# 9.75 hours to run
results_bsa <- mainAnalysis(bsa.data)

########################## Figures - Main Paper ##############################
## To apply par() settings within .eps:
##    1. Set par() settings in postscript embedding procedure
##    2. (ME) Create plots outside of postscript and then call them again inside to save
library(cowplot) # for 'plot_grid'
library(gridGraphics) # for 'plot_grid' functionality with regular plots

## Figure 1 ##
dev.off()
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
# Panel A #
subset_4000 <- subset(gel.data, xor(well == 17, xor(well == 61, well == 92)) & freq == 4000 & tray == 5 & time_index <= 24)
subset_8000 <- subset(gel.data, xor(well == 17, xor(well == 61, well == 92)) & freq == 8000 & tray == 5 & time_index <= 24)
X_a <- matrix(subset_4000[, "res"], 24, 3)
Y_a <- matrix(subset_8000[, "res"], 24, 3)
matplot(X_a, ylab = NA, ylim = c(min(Y_a)-100, max(X_a) + 400), col = c(5, 6, "gold4"), type = "l", lty = 1, lwd = 2)
lines(1:24, Y_a[, 1], col = 5, type = "l", lty = 2, lwd = 2)
lines(1:24, Y_a[, 2], col = 6, type = "l", lty = 2, lwd = 2)
lines(1:24, Y_a[, 3], col = "gold4", type = "l", lty = 2, lwd = 2)
mtext("Resistance (ohms)", side = 2, line = 2.3)
mtext("Time Index (~20 min.)", side = 1, line = 2.2)
legend("topright", legend = c(unique(subset_4000$cell), "@ 4000 Hz", "@ 8000 Hz"), col = c(5, 6, "gold4", 1, 1), lty = c(1,1,1,1, 2), lwd = 2, ncol = 2, bty = 'n', cex = .7)
m1_1 <- recordPlot()
# Panel B #
d <- subset(gel.data, well == 17 & tray == 5 & time_index <= 64 & freq >= 250)
X_b <- matrix(d[, 11], 64, 9)
matplot(X_b, xlab = "Time Index (~20 min.)", ylab = "Resistance (ohms)", ylim = c(0, max(X_b)+4000), col = c(1, 2, 3, 4, 5, 6, 7, 8, "orange"), type = "l", lty = 1, lwd = 2)
legend("topleft", legend = paste(unique(d$freq), "Hz", sep = ""), col = c(1, 2, 3, 4, 5, 6, 7, 8, "orange"), lty = 1, lwd = 2, ncol = 2, bty = "n", cex = .7)
m2_1 <- recordPlot()
# Put together and save # 
setEPS()
postscript("figure1.eps", width = 10, height = 4, family = "Times")
plot_grid(m1_1, m2_1, labels = c("(A)", "(B)"), label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

## Figure 2 ##
# Picture pulled from internet: http://biophysics.com/ecismodel.php #
# Conversion to .eps performed externally #

## Figure 3 ## 
dev.off()
fig3_data <- results_gel$pairfeat$best_feat_data[,,1] # data associated with best features for classification tree method
fig3_lab <- createLab(results_gel$pairfeat$best_feats[1, 1])
fig3_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
fig3_a_parms <- getExample(results_gel$cells, fig3_data, class_type = "tree", train = TRUE)
plot(fig3_data[, 1], fig3_data[, 2], col = rep(fig3_cols, each = 14),
     xlab = fig3_lab[1], ylab = fig3_lab[2], pch = fig3_a_parms$pch_vals)
contour(x = fig3_a_parms$nd.x, y = fig3_a_parms$nd.y, z = matrix(fig3_a_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
m1_3 <- recordPlot()
# Panel B #
fig3_b_parms <- getExample(results_gel$cells, fig3_data, class_type = "tree", train = FALSE)
plot(fig3_data[, 1], fig3_data[, 2], col = rep(fig3_cols, each = 14),
     xlab = fig3_lab[1], ylab = fig3_lab[2], pch = fig3_b_parms$pch_vals)
contour(x = fig3_b_parms$nd.x, y = fig3_b_parms$nd.y, z = matrix(fig3_b_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
m2_3 <- recordPlot()
# Legend #
dev.off()
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    omi = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0)
)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend = unique(results_gel$cells),
       fill = fig3_cols, border = fig3_cols, bty = "n")
m3_3 <- recordPlot()
# Put together and save #
setEPS()
postscript("figure3.eps", width = 12, height = 4, family = "Times")
plot_grid(m1_3, m2_3, m3_3, labels = c("(A)", "(B)", NULL), ncol = 3, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8, rel_widths = c(5, 5, 1))
dev.off()

## Figure 4 ## 
dev.off()
fig4_data <- results_gel$pair_rdascan$best_feat_data
fig4_lab <- createLab(results_gel$pair_rdascan$best_feats[1])
fig4_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
fig4_a_parms <- getExample(results_gel$cells, fig4_data, class_type = "rda", rho_1 = results_gel$pair_rdascan$best_feats[4], train = TRUE)
plot(fig4_data[, 1], fig4_data[, 2], col = rep(fig4_cols, each = 14),
     xlab = fig4_lab[1], ylab = fig4_lab[2], pch = fig4_a_parms$pch_vals)
contour(x = fig4_a_parms$nd.x, y = fig4_a_parms$nd.y, z = matrix(fig4_a_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
m1_4 <- recordPlot()
# Panel B #
fig4_b_parms <- getExample(results_gel$cells, fig4_data, class_type = "rda", rho_1 = results_gel$pair_rdascan$best_feats[4], train = FALSE)
plot(fig4_data[, 1], fig4_data[, 2], col = rep(fig4_cols, each = 14),
     xlab = fig4_lab[1], ylab = fig4_lab[2], pch = fig4_b_parms$pch_vals)
contour(x = fig4_b_parms$nd.x, y = fig4_b_parms$nd.y, z = matrix(fig4_b_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
m2_4 <- recordPlot()
# Legend #
dev.off()
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    omi = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 0)
)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend = unique(results_gel$cells),
       fill = fig4_cols, border = fig4_cols, bty = "n")
m3_4 <- recordPlot()
# Put together and save #
setEPS()
postscript("figure4.eps", width = 12, height = 4, family = "Times")
plot_grid(m1_4, m2_4, m3_4, labels = c("(A)", "(B)", NULL), ncol = 3, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8, rel_widths = c(5, 5, 1))
dev.off()

########################## Tables - Main Paper ##############################
library(xtable)
bold <- function(x){
  paste0('{\\bfseries ', x, '}') }

## Table 1 ##
feature_names_1 <- c(createLab(results_gel$singlefeat$best_feats[1, 1]), paste(createLab(results_gel$pairfeat$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(results_gel$pairfeat$best_feats[1, 1])[2]),
                     paste(createLab(results_gel$triofeat$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[1, 1])[2], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[1, 1])[3]))
accuracy_1 <- c(results_gel$singlefeat$best_feats[1, 2], results_gel$pairfeat$best_feats[1, 2], results_gel$triofeat$best_feats[1, 2])
error_1 <- c(results_gel$singlefeat$best_feats[1, 3], results_gel$pairfeat$best_feats[1, 3], results_gel$triofeat$best_feats[1, 3])
table1 <- data.frame(as.character(c(1, 2, 3)), feature_names_1, accuracy_1, error_1)
colnames(table1) <- c("Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table1, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{6cm}", "C{2.5cm}", "C{2.5cm}")), # first align character is alignment of overall table
      sanitize.colnames.function = bold,  include.rownames = FALSE, hline.after = c(0, 1, 2), sanitize.text = force, size = "small", file = "table1.tex", floating = FALSE)


## Table 2 ##
feature_names_2 <- c(paste(createLab(results_gel$pair_rdascan$best_feats[1])[1], "\\hspace{3cm}", createLab(results_gel$pair_rdascan$best_feats[1])[2]),
                     paste(createLab(results_gel$trio_rdascan$best_feats[1])[1], "\\hspace{3cm}", createLab(results_gel$trio_rdascan$best_feats[1])[2], "\\hspace{3cm}", createLab(results_gel$trio_rdascan$best_feats[1])[3]))
accuracy_2 <- as.numeric(c(results_gel$pair_rdascan$best_feats[2], results_gel$trio_rdascan$best_feats[2]))
error_2 <- as.numeric(c(results_gel$pair_rdascan$best_feats[3], results_gel$trio_rdascan$best_feats[3]))
rho_1s <- c(results_gel$pair_rdascan$best_feats[4], results_gel$trio_rdascan$best_feats[4])
table2 <- data.frame(as.character(c(2, 3)), rho_1s, feature_names_2, accuracy_2, error_2)
colnames(table2) <- c("Feature Space Dimension", "$\\bm{\\rho_1}$", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table2, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{1cm}", "C{5.8cm}", "C{2.5cm}", "C{2.5cm}")), sanitize.colnames.function = bold,  include.rownames = FALSE, # first align character is alignment of overall table
      hline.after = c(0, 1), size = "small", file = "table2.tex", sanitize.text = force, floating = FALSE)

########################## Tables - Supplemental Materials ##############################
library(xtable)
bold <- function(x){
  paste0('{\\bfseries ', x, '}') }

## Supplemental Table 1 ##
feature_names_1_s <- c(createLab(results_gel$singlefeat$best_feats[2, 1]), paste(createLab(results_gel$pairfeat$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(results_gel$pairfeat$best_feats[2, 1])[2]),
                     paste(createLab(results_gel$triofeat$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[2, 1])[2], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[2, 1])[3]))
accuracy_1_s <- c(results_gel$singlefeat$best_feats[2, 4], results_gel$pairfeat$best_feats[2, 4], results_gel$triofeat$best_feats[2, 4])
error_1_s <- c(results_gel$singlefeat$best_feats[2, 5], results_gel$pairfeat$best_feats[2, 5], results_gel$triofeat$best_feats[2, 5])
table1_s <- data.frame(as.character(c(1, 2, 3)), feature_names_1_s, accuracy_1_s, error_1_s)
colnames(table1_s) <- c("Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table1_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{6cm}", "C{2.5cm}", "C{2.5cm}")), # first align character is alignment of overall table
      sanitize.colnames.function = bold,  include.rownames = FALSE, sanitize.text = force, hline.after = c(0, 1, 2), size = "small", file = "supplementaltable1.tex", floating = FALSE)

## Supplemental Table 2 ##
feature_names_2_s <- c(createLab(results_gel$singlefeat$best_feats[3, 1]), paste(createLab(results_gel$pairfeat$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(results_gel$pairfeat$best_feats[3, 1])[2]),
                       paste(createLab(results_gel$triofeat$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[3, 1])[2], "\\hspace{3cm}", createLab(results_gel$triofeat$best_feats[3, 1])[3]))
accuracy_2_s <- c(results_gel$singlefeat$best_feats[3, 6], results_gel$pairfeat$best_feats[3, 6], results_gel$triofeat$best_feats[3, 6])
error_2_s <- c(results_gel$singlefeat$best_feats[3, 7], results_gel$pairfeat$best_feats[3, 7], results_gel$triofeat$best_feats[3, 7])
table2_s <- data.frame(as.character(c(1, 2, 3)), feature_names_2_s, accuracy_2_s, error_2_s)
colnames(table2_s) <- c("Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table2_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{6cm}", "C{2.5cm}", "C{2.5cm}")), # first align character is alignment of overall table
      sanitize.colnames.function = bold,  include.rownames = FALSE, sanitize.text = force, hline.after = c(0, 1, 2), size = "small", file = "supplementaltable2.tex", floating = FALSE)


## Supplemental Table 3 ##
methods_3_s <- c("Tree", "", "", paste("RDA $(\\rho_1=",results_bsa$pair_rdascan$best_feats[4],")$"), paste("$(\\rho_1=",results_bsa$trio_rdascan$best_feats[4],")$"),
                 "LDA", "", "", "QDA", "", "")
feature_names_3_s <- c(createLab(results_bsa$singlefeat$best_feats[1, 1]), paste(createLab(results_bsa$pairfeat$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(results_bsa$pairfeat$best_feats[1, 1])[2]),
                       paste(createLab(results_bsa$triofeat$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[1, 1])[2], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[1, 1])[3]),
                       paste(createLab(results_bsa$pair_rdascan$best_feats[1])[1], "\\hspace{3cm}", createLab(results_bsa$pair_rdascan$best_feats[1])[2]),
                       paste(createLab(results_bsa$trio_rdascan$best_feats[1])[1], "\\hspace{3cm}", createLab(results_bsa$trio_rdascan$best_feats[1])[2], "\\hspace{3cm}", createLab(results_bsa$trio_rdascan$best_feats[1])[3]),
                      createLab(results_bsa$singlefeat$best_feats[2, 1]), paste(createLab(results_bsa$pairfeat$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(results_bsa$pairfeat$best_feats[2, 1])[2]),
                      paste(createLab(results_bsa$triofeat$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[2, 1])[2], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[2, 1])[3]),
                      createLab(results_bsa$singlefeat$best_feats[3, 1]), paste(createLab(results_bsa$pairfeat$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(results_bsa$pairfeat$best_feats[3, 1])[2]),
                      paste(createLab(results_bsa$triofeat$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[3, 1])[2], "\\hspace{3cm}", createLab(results_bsa$triofeat$best_feats[3, 1])[3]))
accuracy_3_s <- as.numeric(c(results_bsa$singlefeat$best_feats[1, 2], results_bsa$pairfeat$best_feats[1, 2], results_bsa$triofeat$best_feats[1, 2],
                  as.numeric(c(results_bsa$pair_rdascan$best_feats[2], results_bsa$trio_rdascan$best_feats[2])),
                  results_bsa$singlefeat$best_feats[2, 4], results_bsa$pairfeat$best_feats[2, 4], results_bsa$triofeat$best_feats[2, 4],
                  results_bsa$singlefeat$best_feats[3, 6], results_bsa$pairfeat$best_feats[3, 6], results_bsa$triofeat$best_feats[3, 6]))
error_3_s <- as.numeric(c(results_bsa$singlefeat$best_feats[1, 3], results_bsa$pairfeat$best_feats[1, 3], results_bsa$triofeat$best_feats[1, 3],
               as.numeric(c(results_bsa$pair_rdascan$best_feats[3], results_bsa$trio_rdascan$best_feats[3])),
               results_bsa$singlefeat$best_feats[2, 5], results_bsa$pairfeat$best_feats[2, 5], results_bsa$triofeat$best_feats[2, 5],
               results_bsa$singlefeat$best_feats[3, 7], results_bsa$pairfeat$best_feats[3, 7], results_bsa$triofeat$best_feats[3, 7]))
table3_s <- data.frame(methods_3_s, as.character(c(1, 2, 3, 2, 3, rep(c(1,2,3), 2))), feature_names_3_s, accuracy_3_s, error_3_s)
colnames(table3_s) <- c("Method", "Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table3_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{2cm}", "C{5.8cm}", "C{2.5cm}", "C{2.5cm}")), sanitize.colnames.function = bold,  include.rownames = FALSE, # first align character is alignment of overall table
      hline.after = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), sanitize.text = force, size = "small", file = "supplementaltable3.tex", floating = FALSE, sanitize.text.function=function(x){x})

## Supplemental Table 4 ##
methods_4_s <- c("Tree", "", "", "LDA", "", "", "QDA", "", "")
feature_names_4_s <- c(createLab(single_wtis$best_feats[1, 1]), paste(createLab(pair_wtis$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(pair_wtis$best_feats[1, 1])[2]),
                       paste(createLab(trio_wtis$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(trio_wtis$best_feats[1, 1])[2], "\\hspace{3cm}", createLab(trio_wtis$best_feats[1, 1])[3]),
                       createLab(single_wtis$best_feats[2, 1]), paste(createLab(pair_wtis$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(pair_wtis$best_feats[2, 1])[2]),
                       paste(createLab(trio_wtis$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(trio_wtis$best_feats[2, 1])[2], "\\hspace{3cm}", createLab(trio_wtis$best_feats[2, 1])[3]),
                       createLab(single_wtis$best_feats[3, 1]), paste(createLab(pair_wtis$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(pair_wtis$best_feats[3, 1])[2]),
                       paste(createLab(trio_wtis$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(trio_wtis$best_feats[3, 1])[2], "\\hspace{3cm}", createLab(trio_wtis$best_feats[3, 1])[3]))
accuracy_4_s <- as.numeric(c(single_wtis$best_feats[1, 2], pair_wtis$best_feats[1, 2], trio_wtis$best_feats[1, 2],
                             single_wtis$best_feats[2, 4], pair_wtis$best_feats[2, 4], trio_wtis$best_feats[2, 4],
                             single_wtis$best_feats[3, 6], pair_wtis$best_feats[3, 6], trio_wtis$best_feats[3, 6]))
error_4_s <- as.numeric(c(single_wtis$best_feats[1, 3], pair_wtis$best_feats[1, 3], trio_wtis$best_feats[1, 3],
                          single_wtis$best_feats[2, 5], pair_wtis$best_feats[2, 5], trio_wtis$best_feats[2, 5],
                          single_wtis$best_feats[3, 7], pair_wtis$best_feats[3, 7], trio_wtis$best_feats[3, 7]))
table4_s <- data.frame(methods_4_s, as.character(rep(c(1,2,3), 3)), feature_names_4_s, accuracy_4_s, error_4_s)
colnames(table4_s) <- c("Method", "Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table4_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{2cm}", "C{5.8cm}", "C{2.5cm}", "C{2.5cm}")), sanitize.colnames.function = bold,  include.rownames = FALSE, # first align character is alignment of overall table
      hline.after = c(0, 1, 2, 3, 4, 5, 6, 7, 8), sanitize.text = force, size = "small", file = "supplementaltable4.tex", floating = FALSE)

## Supplemental Table 5 ##
methods_5_s <- c("Tree", "", "", "LDA", "", "", "QDA", "", "")
feature_names_5_s <- c(createLab(single_wtis_just_tis$best_feats[1, 1]), paste(createLab(pair_wtis_just_tis$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(pair_wtis_just_tis$best_feats[1, 1])[2]),
                       paste(createLab(trio_wtis_just_tis$best_feats[1, 1])[1], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[1, 1])[2], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[1, 1])[3]),
                       createLab(single_wtis_just_tis$best_feats[2, 1]), paste(createLab(pair_wtis_just_tis$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(pair_wtis_just_tis$best_feats[2, 1])[2]),
                       paste(createLab(trio_wtis_just_tis$best_feats[2, 1])[1], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[2, 1])[2], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[2, 1])[3]),
                       createLab(single_wtis_just_tis$best_feats[3, 1]), paste(createLab(pair_wtis_just_tis$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(pair_wtis_just_tis$best_feats[3, 1])[2]),
                       paste(createLab(trio_wtis_just_tis$best_feats[3, 1])[1], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[3, 1])[2], "\\hspace{3cm}", createLab(trio_wtis_just_tis$best_feats[3, 1])[3]))
accuracy_5_s <- as.numeric(c(single_wtis_just_tis$best_feats[1, 2], pair_wtis_just_tis$best_feats[1, 2], trio_wtis_just_tis$best_feats[1, 2],
                             single_wtis_just_tis$best_feats[2, 4], pair_wtis_just_tis$best_feats[2, 4], trio_wtis_just_tis$best_feats[2, 4],
                             single_wtis_just_tis$best_feats[3, 6], pair_wtis_just_tis$best_feats[3, 6], trio_wtis_just_tis$best_feats[3, 6]))
error_5_s <- as.numeric(c(single_wtis_just_tis$best_feats[1, 3], pair_wtis_just_tis$best_feats[1, 3], trio_wtis_just_tis$best_feats[1, 3],
                          single_wtis_just_tis$best_feats[2, 5], pair_wtis_just_tis$best_feats[2, 5], trio_wtis_just_tis$best_feats[2, 5],
                          single_wtis_just_tis$best_feats[3, 7], pair_wtis_just_tis$best_feats[3, 7], trio_wtis_just_tis$best_feats[3, 7]))
table5_s <- data.frame(methods_5_s, as.character(rep(c(1,2,3), 3)), feature_names_5_s, accuracy_5_s, error_5_s)
colnames(table5_s) <- c("Method", "Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table5_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{2cm}", "C{5.8cm}", "C{2.5cm}", "C{2.5cm}")), sanitize.colnames.function = bold,  include.rownames = FALSE, # first align character is alignment of overall table
      hline.after = c(0, 1, 2, 3, 4, 5, 6, 7, 8), sanitize.text = force, size = "small", file = "supplementaltable5.tex", floating = FALSE)

## Supplemental Table 6 ##
feature_names_6_s <- c(createLab(single_knn[1, 1]), paste(createLab(pair_knn[1, 1])[1], "\\hspace{3cm}", createLab(pair_knn[1, 1])[2]),
                     paste(createLab(trio_knn[1, 1])[1], "\\hspace{3cm}", createLab(trio_knn[1, 1])[2], "\\hspace{3cm}", createLab(trio_knn[1, 1])[3]))
accuracy_6_s <- c(single_knn[1, 2], pair_knn[1, 2], trio_knn[1, 2])
error_6_s <- c(single_knn[1, 3], pair_knn[1, 3], trio_knn[1, 3])
table6_s <- data.frame(as.character(c(1, 2, 3)), feature_names_6_s, accuracy_6_s, error_6_s)
colnames(table6_s) <- c("Feature Space Dimension", "Selected Classification Feature(s)", "Out-of-Sample Classification Accuracy", "Approximate Standard Error")
print(xtable(table6_s, type = "latex", digits = 3, align = c("c", "C{2cm}", "C{6cm}", "C{2.5cm}", "C{2.5cm}")), # first align character is alignment of overall table
      sanitize.colnames.function = bold,  include.rownames = FALSE, hline.after = c(0, 1, 2), sanitize.text = force, size = "small", file = "supplementaltable6.tex", floating = FALSE)


################### Figures - Supplemental Materials #####################
## To apply par() settings within .eps:
##    1. Set par() settings in postscript embedding procedure
##    2. (ME) Create plots outside of postscript and then call them again inside to save
library(cowplot) # for 'plot_grid'
library(gridGraphics) # for 'plot_grid' functionality with regular plots
library(whitening) # for 'whiten' to decorrelate variables

## Supplemental Figure 1 ##
dev.off()
rho_vals <- c(0, results_gel$rhos, 1) # include QDA (rho_1 = 0) and LDA (rho_1 = 1) results as well
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
max_pair_rda <- c(results_gel$pairfeat$best_feats[3, 6], apply(results_gel$pair_rdascan$class_rates_RDA, 3, function(x) max(x[, 2])), results_gel$pairfeat$best_feats[2, 4])
plot(rho_vals, max_pair_rda, xlab = expression(rho[1]), ylab = "Out-of-Sample Classification Accuracy", 
     main = "Two Dimensions", type = "o", col = 1, lwd = 2, ylim = c(as.numeric(min(.9, min(max_pair_rda))), 1))
points(rho_vals[which.max(max_pair_rda)], max(max_pair_rda), col = "green", pch = 17, cex = 2)
m1_1s <- recordPlot()
# Panel B #
max_trio_rda <- c(results_gel$triofeat$best_feats[3, 6], apply(results_gel$trio_rdascan$class_rates_RDA, 3, function(x) max(x[, 2])), results_gel$triofeat$best_feats[2, 4])
plot(rho_vals, max_trio_rda, xlab = expression(rho[1]), ylab = "Out-of-Sample Classification Accuracy",
     main = "Three Dimensions", type = "o", col = 1, lwd = 2, ylim = c(as.numeric(min(.9, min(max_trio_rda))), 1))
points(rho_vals[which.max(max_trio_rda)], max(max_trio_rda), col = "green", pch = 17, cex = 2)
m2_1s <- recordPlot()
# Put together and save #
setEPS()
postscript("supplementalfigure1.eps", width = 12, height = 4, family = "Times")
plot_grid(m1_1s, m2_1s, labels = c("(A)", "(B)"), ncol = 2, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

## Supplemental Figure 2 ##
dev.off()
fig2_s_data <- results_gel$pair_rdascan$best_feat_data # data associated with best features for RDA method
fig2_s_lab <- createLab(results_gel$pair_rdascan$best_feats[1])
fig2_s_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0), # move axis labels closer to axis
    mar = par()$mar + c(0,0,0,3)
)
fig2_s_parms <- list(getExample(results_gel$cells, fig2_s_data, class_type = "rda", rho_1 = 0.1, train = FALSE),
                     getExample(results_gel$cells, fig2_s_data, class_type = "rda", rho_1 = 0.3, train = FALSE),
                     getExample(results_gel$cells, fig2_s_data, class_type = "rda", rho_1 = 0.7, train = FALSE))
plot(fig2_s_data[, 1], fig2_s_data[, 2], col = rep(fig2_s_cols, each = 14),
     xlab = fig2_s_lab[1], ylab = fig2_s_lab[2], pch = fig2_s_parms[[1]]$pch_vals)
contour(x = fig2_s_parms[[1]]$nd.x, y = fig2_s_parms[[1]]$nd.y, z = matrix(fig2_s_parms[[1]]$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
contour(x = fig2_s_parms[[2]]$nd.x, y = fig2_s_parms[[2]]$nd.y, z = matrix(fig2_s_parms[[2]]$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE, lty = 2)
contour(x = fig2_s_parms[[3]]$nd.x, y = fig2_s_parms[[3]]$nd.y, z = matrix(fig2_s_parms[[3]]$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE, lty = 3)
legend(max(fig2_s_data[, 1]) + 100, max(fig2_s_data[, 2]) - 2800, legend = unique(results_gel$cells), fill = fig2_s_cols, border = fig2_s_cols, cex = .65, ncol = 1, bty = "n")
legend(max(fig2_s_data[, 1]) + 100, max(fig2_s_data[, 2]) + 1000, legend = c(0.1, 0.3, 0.7), lty = c(2, 1, 3), ncol = 1, cex = .7, title = expression(paste(rho[1])), bty = "n")
m1_2s <- recordPlot()
setEPS()
postscript("supplementalfigure2.eps", width = 7.5, height = 6, family = "Times")
m1_2s
dev.off()

## Supplemental Figure 3 ##
dev.off()
fig3_s_data <- results_gel$pairfeat$best_feat_data[,,2] # data associated with best features for LDA method
fig3_s_lab <- createLab(results_gel$pairfeat$best_feats[2, 1])
fig3_s_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0), # move axis labels closer to axis
    mar = par()$mar + c(0, 0, 0, 5) # make room for legends
)
fig3_s_a_parms <- getExample(results_gel$cells, fig3_s_data, class_type = "lda", train = TRUE)
plot(fig3_s_data[, 1], fig3_s_data[, 2], col = rep(fig3_s_cols, each = 14),
     xlab = fig3_s_lab[1], ylab = fig3_s_lab[2], pch = fig3_s_a_parms$pch_vals)
contour(x = fig3_s_a_parms$nd.x, y = fig3_s_a_parms$nd.y, z = matrix(fig3_s_a_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig3_s_data[, 1]) + 200, y = max(fig3_s_data[, 2]) + 800, legend = unique(results_gel$cells),
       fill = fig3_s_cols, border = fig3_s_cols, bty = "n", cex = 0.8)
m1_3s <- recordPlot()
# Panel B #
fig3_s_b_parms <- getExample(results_gel$cells, fig3_s_data, class_type = "lda", train = FALSE)
plot(fig3_s_data[, 1], fig3_s_data[, 2], col = rep(fig3_s_cols, each = 14),
     xlab = fig3_s_lab[1], ylab = fig3_s_lab[2], pch = fig3_s_b_parms$pch_vals)
contour(x = fig3_s_b_parms$nd.x, y = fig3_s_b_parms$nd.y, z = matrix(fig3_s_b_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig3_s_data[, 1]) + 200, y = max(fig3_s_data[, 2]) + 800, legend = unique(results_gel$cells),
       fill = fig3_s_cols, border = fig3_s_cols, bty = "n", cex = 0.8)
m2_3s <- recordPlot()
# Put together and save #
setEPS()
postscript("supplementalfigure3.eps", width = 6, height = 7, family = "Times")
plot_grid(m1_3s, m2_3s, labels = c("(A)", "(B)"), nrow = 2, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

## Supplemental Figure 4 ##
dev.off()
#fig4_s_data <- results_gel$pairfeat$best_feat_data[,,3] # data associated with best features for QDA method
fig4_s_data <- whiten(results_gel$pairfeat$best_feat_data[,,3], center = TRUE, method = "ZCA")
fig4_s_lab <- createLab(results_gel$pairfeat$best_feats[3, 1])
fig4_s_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0), # move axis labels closer to axis
    mar = par()$mar + c(0, 0, 0, 5) # make room for legends
)
fig4_s_a_parms <- getExample(results_gel$cells, fig4_s_data, class_type = "qda", train = TRUE)
plot(fig4_s_data[, 1], fig4_s_data[, 2], col = rep(fig4_s_cols, each = 14),
     xlab = fig4_s_lab[1], ylab = fig4_s_lab[2], pch = fig4_s_a_parms$pch_vals)
contour(x = fig4_s_a_parms$nd.x, y = fig4_s_a_parms$nd.y, z = matrix(fig4_s_a_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig4_s_data[, 1]) + 0.3, y = max(fig4_s_data[, 2]) + 0.3, legend = unique(results_gel$cells),
       fill = fig4_s_cols, border = fig4_s_cols, bty = "n", cex = 0.8)
m1_4s <- recordPlot()
# Panel B #
fig4_s_b_parms <- getExample(results_gel$cells, fig4_s_data, class_type = "qda", train = FALSE)
plot(fig4_s_data[, 1], fig4_s_data[, 2], col = rep(fig4_s_cols, each = 14),
     xlab = fig4_s_lab[1], ylab = fig4_s_lab[2], pch = fig4_s_b_parms$pch_vals)
contour(x = fig4_s_b_parms$nd.x, y = fig4_s_b_parms$nd.y, z = matrix(fig4_s_b_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig4_s_data[, 1]) + 0.3, y = max(fig4_s_data[, 2]) + 0.3, legend = unique(results_gel$cells),
       fill = fig4_s_cols, border = fig4_s_cols, bty = "n", cex = 0.8)
m2_4s <- recordPlot()
# Put together and save #
setEPS()
postscript("supplementalfigure4.eps", width = 6, height = 7, family = "Times")
plot_grid(m1_4s, m2_4s, labels = c("(A)", "(B)"), nrow = 2, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

## Supplemental Figure 5 ##
dev.off()
rho_vals_bsa <- c(0, results_bsa$rhos, 1) # include QDA (rho_1 = 0) and LDA (rho_1 = 1) results as well
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
max_pair_rda_bsa <- c(results_bsa$pairfeat$best_feats[3, 6], apply(results_bsa$pair_rdascan$class_rates_RDA, 3, function(x) max(x[, 2])), results_bsa$pairfeat$best_feats[2, 4])
plot(rho_vals_bsa, max_pair_rda_bsa, xlab = expression(rho[1]), ylab = "Out-of-Sample Classification Accuracy", 
     main = "Two Dimensions", type = "o", col = 1, lwd = 2, ylim = c(as.numeric(min(.9, min(max_pair_rda_bsa))), 1))
points(rho_vals_bsa[which.max(max_pair_rda_bsa)], max(max_pair_rda_bsa), col = "green", pch = 17, cex = 2)
m1_5s <- recordPlot()
# Panel B #
max_trio_rda_bsa <- c(results_bsa$triofeat$best_feats[3, 6], apply(results_bsa$trio_rdascan$class_rates_RDA, 3, function(x) max(x[, 2])), results_bsa$triofeat$best_feats[2, 4])
plot(rho_vals_bsa, max_trio_rda_bsa, xlab = expression(rho[1]), ylab = "Out-of-Sample Classification Accuracy", 
     main = "Three Dimensions", type = "o", col = 1, lwd = 2, ylim = c(as.numeric(min(.9, min(max_trio_rda_bsa))), 1))
points(rho_vals_bsa[which.max(max_trio_rda_bsa)], max(max_trio_rda_bsa), col = "green", pch = 17, cex = 2)
m2_5s <- recordPlot()
# Put together and save #
setEPS()
postscript("supplementalfigure5.eps", width = 12, height = 4, family = "Times")
plot_grid(m1_5s, m2_5s, labels = c("(A)", "(B)"), ncol = 2, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

## Supplemental Figure 6 ##
dev.off()
#fig6_s_data <- results_bsa$pair_rdascan$best_feat_data # data associated with best features for RDA method
fig6_s_data <- whiten(results_bsa$pair_rdascan$best_feat_data, center = TRUE, method = "ZCA")
fig6_s_lab <- createLab(results_bsa$pair_rdascan$best_feats[1])
fig6_s_cols <- c(1, "slategray1", 3, 4, 5, 6, 7, 8, "orange", "pink", "darkorange3", "forestgreen", "gold4", "blueviolet", "firebrick")
# Panel A #
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(1, 1, 0, 0), # move plot to the right and up
    mar = par()$mar + c(0, 0, 0, 5), # make room for legends
    mgp = c(2.5, 1, 0) # move axis labels closer to axis
)
fig6_s_a_parms <- getExample(results_bsa$cells, fig6_s_data, class_type = "rda", rho_1 = as.numeric(results_bsa$pair_rdascan$best_feats[4]), train = TRUE)
plot(fig6_s_data[, 1], fig6_s_data[, 2], col = rep(fig6_s_cols, each = 14),
     xlab = fig6_s_lab[1], ylab = fig6_s_lab[2], pch = fig6_s_a_parms$pch_vals)
contour(x = fig6_s_a_parms$nd.x, y = fig6_s_a_parms$nd.y, z = matrix(fig6_s_a_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig6_s_data[, 1]) + 0.3, y = max(fig6_s_data[, 2]) + 0.3, legend = unique(results_bsa$cells),
       fill = fig6_s_cols, border = fig6_s_cols, bty = "n", cex = 0.8)
m1_6s <- recordPlot()
# Panel B #
fig6_s_b_parms <- getExample(results_bsa$cells, fig6_s_data, class_type = "rda", rho_1 = as.numeric(results_bsa$pair_rdascan$best_feats[4]), train = FALSE)
plot(fig6_s_data[, 1], fig6_s_data[, 2], col = rep(fig6_s_cols, each = 14),
     xlab = fig6_s_lab[1], ylab = fig6_s_lab[2], pch = fig6_s_b_parms$pch_vals)
contour(x = fig6_s_b_parms$nd.x, y = fig6_s_b_parms$nd.y, z = matrix(fig6_s_b_parms$prd, nrow = 300, ncol = 300), 
        levels = c(1:15), add = TRUE, drawlabels = FALSE)
legend(x = max(fig6_s_data[, 1]) + 0.3, y = max(fig6_s_data[, 2]) + 0.3, legend = unique(results_bsa$cells),
       fill = fig6_s_cols, border = fig6_s_cols, bty = "n", cex = 0.8)
m2_6s <- recordPlot()
# Put together and save #
setEPS()
postscript("supplementalfigure6.eps", width = 6, height = 7, family = "Times")
plot_grid(m1_6s, m2_6s, labels = c("(A)", "(B)"), nrow = 2, label_x = 0, hjust = -0.2, label_y = 0, vjust = -0.8)
dev.off()

