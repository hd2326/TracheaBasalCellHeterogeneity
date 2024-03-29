# Raw umi matrix is avaliable at GEO under accession number: GSE134064
# tpm matrix is generated using: tpm <- log2(t(t(umi)/(colSums(umi)/1e6)) + 1)

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))

pdf("Fig1C.pdf", width = 10, height = 6)
layout(matrix(1:2, 1, 2), widths = c(6, 4))
par(mar = c(5, 5, 1, 1))
colors <- c(rgb(43, 120, 182, maxColorValue = 255),
            rgb(247, 126, 11, maxColorValue = 255),
            rgb(59, 160, 30, maxColorValue = 255),
            rgb(206, 37, 41, maxColorValue = 255),
            rgb(143, 103, 191, maxColorValue = 255),
            rgb(134, 86, 74, maxColorValue = 255))
plot(tsne[clusters %in% 1:6, ], pch = 16, col = colors[clusters[clusters %in% 1:6]], xlab = "tSNE-1", ylab = "tSNE-2", cex.axis = 2, cex.lab = 2, main = "")
par(mar = c(0, 0, 0, 0))
plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = "", ylab = "", main = "")
legend(0, 0.9, legend = c("Basal-1, n=1939",
                             "Basal-2, n=1669",
                             "Secretory, n=243",
                             "Squamous, n=139",
                             "Mes-like, n=96",
                             "Proliferating, n=90"), y.intersp = 2, fill = colors, text.font = 2, cex = 1.5, bty = "n", border = NA)
dev.off()
#Figure 1C

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")

pdf("Fig1D.pdf", width = 12, height = 16)
par(mfrow = c(4, 3), mar = c(0, 0, 5, 0))
for (x in c("Epcam", "Nkx2-1", "Mki67",
            "Trp63", "Krt5", "Krt15", 
            "Krt14", "Krt13", "Krt4",
            "Scgb1a1", "Muc5b", "Vim")){
  colors <- tpm[x, clusters %in% 1:6]
  colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 4)}
dev.off()
#Figure 1D

clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
e <- tpm[c("Mki67", "Ccna2", "Aurkb", "Ube2c", "Pbk", 
           "Nusap1", "Cd3g", "Gng2", "Ptprcap", "Nkg7", 
           "Trbc1", "Ccl5", "Cd52", "Id2", "Thy1", 
           "Vim", "Tgfbi", "Igfbp2", "Mal", "Ly6d", 
           "Sprr2a3", "Serpinb2", "Lypd3", "Dmkn", "Krt13", 
           "Krt6a", "Krt14", "Tff2", "Aldh1a1", "Creb3l1", 
           "Pon1", "Lyz2", "Bpifb1", "Muc5b", "Scgb3a1", 
           "Scgb3a2", "Scgb1a1", "Pdpn", "Aqp3", "Aqp5", 
           "Trp63", "Krt5", "Nkx2-1", "Epcam"),
         c(which(clusters == 1), which(clusters == 2),
           which(clusters == 3), which(clusters == 4),
           which(clusters == 5), which(clusters == 6))]

pdf("SFig1C.pdf", width = 14, height = 6)
layout(matrix(1:2, 1, 2), widths = c(2, 12))
par(mar = c(17, 1, 10, 1))
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Grey", "Red"))(20), axes = F, main = "log2(tpm+1)")
axis(side = 1, at = c(0, 1), labels = c(0, floor(max(e)+1)), tick = F, cex.axis = 1)
par(mar = c(4, 4, 4, 1))
image(t(e), col = colorRampPalette(c("Grey", "Red"))(20), axes = F,
      breaks = seq(0, floor(max(e)+1), length.out = 21), main = "")
abline(NULL, NULL, NULL, cumsum(c(0, sum(clusters == 1), sum(clusters == 2),
                                  sum(clusters == 3), sum(clusters == 4),
                                  sum(clusters == 5), sum(clusters == 6)))/ncol(e))
axis(side = 2, line = 0, at = seq(0, 1, length.out = nrow(e)), labels = rownames(e), tick = F, las = 2, cex.axis = 0.5, font.axis = 3)
axis(side = 1, line = 0, at = cumsum(c(sum(clusters == 1), sum(clusters == 2),
                                       sum(clusters == 3), sum(clusters == 4),
                                       sum(clusters == 5), sum(clusters == 6)))/ncol(e), labels = NA, tick = F)
dev.off()
#SFigure 1C

library(vioplot)
tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
colors <- c(rgb(43, 120, 182, maxColorValue = 255),
            rgb(247, 126, 11, maxColorValue = 255),
            rgb(59, 160, 30, maxColorValue = 255),
            rgb(206, 37, 41, maxColorValue = 255),
            rgb(143, 103, 191, maxColorValue = 255),
            rgb(134, 86, 74, maxColorValue = 255))

pdf("SFig2.pdf", width = 15, height = 20)
layout(matrix(c(1, 1, 2, 2, 3, 3,
                4, 4, 5, 5, 6, 6,
                7, 7, 8, 8, 9, 9,
                10, 10, 11, 11, 12, 12,
                13, 13, 14, 14, 15, 15,
                16, 16, 16, 17, 17, 17), 6, 6, byrow = T), heights = c(2, 2, 2, 2, 2, 3))
par(mar = c(4, 5, 4, 1))
for (gene in c("Epcam", "Nkx2-1", "Trp63", "Krt5", "Krt15",
               "Krt14", "Krt13", "Krt4", "Scgb1a1", "Scgb3a2",
               "Muc5b", "Vim", "Thy1", "Mki67", "Aurkb")){
  table <- lapply(1:6, function(i, clusters, x){
    x[clusters == i]
  }, clusters=clusters[clusters %in% 1:6], x=tpm[gene, clusters %in% 1:6])
  vioplot(table, col = colors, names = c("BC1", "BC2", "Sec.", "Sq.", "Mes.", "Pro."),
          ylab = "log2(tpm+1)", cex.lab = 2, cex.axis = 2)
  axis(side = 3, at = 3.5, labels = gene, cex.axis = 3, tick = F, font = 4)}

par(mar = c(4, 5, 8, 1))
bc <- c("Trp63", "Krt5", "Snai2", "Egfr", "Krt15", "Jag2",
        "Dll1", "Ngfr", "Dapl1", "Gpr87", "Dlk2", "Bcam")
z <- t(apply(tpm[bc, clusters %in% 1:6], 1, scale))
z <- colSums(z)/sqrt(nrow(z))
table <- lapply(1:6, function(i, clusters, z) z[clusters == i], clusters=clusters[clusters %in% 1:6], z=z)
vioplot(table, col = colors, names = rep("", 6), cex.axis = 3, main = "Basal Score", cex.main = 4)
axis(side = 1, line = 1, at = 1:6, labels = c("BC1", "BC2", "Sec.", "Sq.", "Mes.", "Pro."), tick = F, cex.axis = 3)

par(mar = c(4, 5, 8, 1))
lm <- c("Krt8", "Scgb1a1", "Scgb3a2", "Muc1", "Krt13", "Krt4", "Lypd2",
        "Notch2", "Notch3", "Perp", "Ly6d", "Lypd3", "Ly6e", "Hes1")
z <- t(apply(tpm[lm, clusters %in% 1:6], 1, scale))
z <- colSums(z)/sqrt(nrow(z))
table <- lapply(1:6, function(i, clusters, z) z[clusters == i], clusters=clusters[clusters %in% 1:6], z=z)
vioplot(table, col = colors, names = rep("", 6), cex.axis = 3, main = "Luminal Score", cex.main = 4)
axis(side = 1, line = 1, at = 1:6, labels = c("BC1", "BC2", "Sec.", "Sq.", "Mes.", "Pro."), tick = F, cex.axis = 3)
dev.off()
#SFigure 2
                
tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
d <- rowMeans(tpm[, clusters == 1]) - rowMeans(tpm[, clusters == 2])
z <- (tpm - rowMeans(tpm))/apply(tpm, 1, sd)
z <- z[!is.na(rowSums(z)), ]
z <- z[rev(c(names(sort(d, decreasing = T))[1:50], rev(names(sort(d, decreasing = F))[1:50]))), 
       c(which(clusters == 1), which(clusters == 2))]

pdf("SFig3A.pdf", width = 14, height = 22)
layout(matrix(1:2, 1, 2, byrow = T), widths = c(2, 4))
par(mar = c(52, 1, 52, 6))
image(matrix(1:100, 100, 1, byrow = T), col = colorRampPalette(c("Blue", "White", "Red"))(100), axes = F, main = "")
axis(side = 3, at = 0.5, labels = "Z-Score", tick = F, cex.axis = 2)
axis(side = 1, at = c(0, 0.5, 1), labels = c(-4, 0, 4), tick = F, cex.axis = 2)

par(mar = c(1, 6, 4, 1))
image(t(z), col = colorRampPalette(c("Blue", "Grey", "Red"))(100), axes = F,
      breaks = seq(-4, 4, length.out = 101), main = "Basal-1 vs. Basal-2", cex.main = 4)
abline(NULL, NULL, NULL, cumsum(c(0, sum(clusters == 1), sum(clusters == 2)))/ncol(z))
axis(side = 2, at = seq(0, 1, length.out = nrow(z)), labels = rownames(z), tick = F, las = 2, cex.axis = 1.5, font.axis = 3)
dev.off()
#SFigure 3A

load("./expression/tpm.rda")
z <- (tpm - rowMeans(tpm))/apply(tpm, 1, sd)
z <- z[!is.na(rowSums(z)), ]
z <- z[c("Klf6", "Maff", "Klf4", "Barx2", "Myc", "Foxq1", "Atf3", "Jun", 
         "Baz1a", "Hivep2", "Egr1", "Jund", "Zfp36", "Nr4a1", "Irf6", "Klf13", "Srf", "Tgif1",
         "Bhlhe40", "Stat3", "Mafg", "Junb", "Nfix", "Tef", "Id1", 
         "Six1", "Tsc22d3", "Mecom", "Isl1", "Cited2"), 
       c(which(clusters == 1), which(clusters == 2))]

pdf("SFig3B.pdf", width = 14, height = 10)
layout(matrix(1:2, 1, 2, byrow = T), widths = c(2, 4))
par(mar = c(26, 1, 20, 6))
image(matrix(1:100, 100, 1, byrow = T), col = colorRampPalette(c("Blue", "White", "Red"))(100), axes = F, main = "")
axis(side = 3, at = 0.5, labels = "Z-Score", tick = F, cex.axis = 2)
axis(side = 1, at = c(0, 0.5, 1), labels = c(-4, 0, 4), tick = F, cex.axis = 2)

par(mar = c(4, 4, 6, 1))
image(t(z), col = colorRampPalette(c("Blue", "Grey", "Red"))(100), axes = F,
      breaks = seq(-4, 4, length.out = 101), main = "Basal-1 vs. Basal-2", cex.main = 4)
abline(NULL, NULL, NULL, cumsum(c(0, sum(clusters == 1), sum(clusters == 2)))/ncol(z))
axis(side = 2, at = seq(0, 1, length.out = nrow(z)), labels = rownames(z), tick = F, las = 2, cex.axis = 1.5, font.axis = 3)
dev.off()
#SFigure 3B

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
tpm <- tpm[!duplicated(rownames(tpm)), ]
p <- -log10(p.adjust(unlist(lapply(1:nrow(tpm), function(i, x, y){
  wilcox.test(x[i, ], y[i, ], alternative = "two.sided")$p.value
}, x=tpm[, clusters == 4], y=tpm[, clusters %in% c(1, 2, 3, 5, 6)]))))
d <- rowMeans(tpm[, clusters == 4]) - rowMeans(tpm[, clusters %in% c(1, 2, 3, 5, 6)])
names(p) <- names(d) <- rownames(tpm)

pdf("SFig4A.pdf", width = 5, height = 5)
par(mar = c(5, 5, 5, 1))
genes <- c("Krt6a", "Krt13", "Serpinb2", "Mal", "Krt14", "Krt4", "Dmkn", "Sprr1a")
plot(d, p, pch = 16, cex = 0.5, xlab = "d", ylab = "-log10(p)", xlim = c(-8, 8),
     cex.axis = 2, cex.lab = 2, main = "Squamous vs. Others", cex.main = 2)
points(d[genes], p[genes], col = 2, pch = 16, cex = 1)
text(d[genes], p[genes], labels = genes, col = 2, cex = 1, pos = 2, font = 4)
dev.off()
#SFigure 4A

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
z <- (tpm - rowMeans(tpm))/apply(tpm, 1, sd)
z <- z[!is.na(rowSums(z)), ]
z <- z[c("Trp53inp2", "Cldn4", "Barx2", "Crip1", "Pdzk1ip1", 
         "Serpinb1a", "Mboat1", "Pmm1", "Cd24a", "Hebp2", 
         "Pkp1", "St3gal4", "Tuba1a", "Ceacam1", "Nab1", 
         "Card19", "Dsc2", "Tmem43", "Cda", "Lypd3", 
         "Igfbp3", "Lad1", "Serpinb10", "Serpinb5", "Sprr1a", 
         "Tmprss4", "Cdkn1a", "Ociad2", "Krt17", "Upk3bl", 
         "Nupr1", "Dmkn", "Tppp3", "Krt4", "Tpm2", 
         "Ltf", "Ly6g6c", "Ly6c1", "Ly6d", "Sprr2a3", 
         "Lgals7", "Ecm1", "Calml3", "Krt14", "Plac8", 
         "Mal", "Cwh43", "Serpinb2", "Krt13", "Krt6a"), 
       c(which(clusters == 1), which(clusters == 2),
         which(clusters == 3), which(clusters == 4),
         which(clusters == 5), which(clusters == 6))]

pdf("SFig4C.pdf", width = 18, height = 12)
layout(matrix(1:2, 1, 2, byrow = T), widths = c(2, 4))
par(mar = c(45, 1, 10, 10))
image(matrix(1:100, 100, 1, byrow = T), col = colorRampPalette(c("Blue", "White", "Red"))(100), axes = F, main = "")
axis(side = 3, at = 0.5, labels = "Z-Score", tick = F, cex.axis = 2)
axis(side = 1, at = c(0, 0.5, 1), labels = c(-4, 0, 4), tick = F, cex.axis = 2)

par(mar = c(1, 4, 1, 1))
image(t(z), col = colorRampPalette(c("Blue", "Grey", "Red"))(100),
      axes = F, breaks = seq(-4, 4, length.out = 101), main = "")
abline(NULL, NULL, NULL, cumsum(c(0, sum(clusters == 1), sum(clusters == 2),
                                  sum(clusters == 3), sum(clusters == 4),
                                  sum(clusters == 5), sum(clusters == 6)))/ncol(z))
axis(side = 2, at = seq(0, 1, length.out = nrow(z)), labels = rownames(z), tick = F, las = 2, cex.axis = 1.5, font.axis = 3)
dev.off()
#SFigure 4C

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")
z <- (tpm - rowMeans(tpm))/apply(tpm, 1, sd)
z <- z[!is.na(rowSums(z)), ]
z <- z[c("Tsc22d3", "Trp63", "Srebf2", "Pax9", "Klf4", 
         "Ikzf2", "Id1", "Hmga1", "Bhlhe40", "Barx2"), 
       c(which(clusters == 1), which(clusters == 2),
         which(clusters == 3), which(clusters == 4),
         which(clusters == 5), which(clusters == 6))]

pdf("SFig4D.pdf", width = 18, height = 5)
layout(matrix(1:2, 1, 2, byrow = T), widths = c(2, 4))
par(mar = c(10, 1, 10, 10))
image(matrix(1:100, 100, 1, byrow = T), col = colorRampPalette(c("Blue", "White", "Red"))(100), axes = F, main = "")
axis(side = 3, at = 0.5, labels = "Z-Score", tick = F, cex.axis = 2)
axis(side = 1, at = c(0, 0.5, 1), labels = c(-4, 0, 4), tick = F, cex.axis = 2)

par(mar = c(4, 4, 4, 1))
image(t(z), col = colorRampPalette(c("Blue", "Grey", "Red"))(100), axes = F,
      breaks = seq(-4, 4, length.out = 101), main = "")
abline(NULL, NULL, NULL, cumsum(c(0, sum(clusters == 1), sum(clusters == 2),
                                  sum(clusters == 3), sum(clusters == 4),
                                  sum(clusters == 5), sum(clusters == 6)))/ncol(z))
axis(side = 2, at = seq(0, 1, length.out = nrow(z)), labels = rownames(z), tick = F, las = 2, cex.axis = 1.5, font.axis = 3)
dev.off()
#SFigure 4D

tsne <- as.matrix(read.csv("./tsne/2_components/projection.csv", row.names=1))
clusters <- as.matrix(read.csv("./clustering/kmeans_9_clusters/clusters.csv", row.names=1))
clusters <- structure(clusters[, 1], names = rownames(clusters))
load("./expression/tpm.rda")

pdf("SFig5.pdf", width = 24, height = 22)
par(mfrow = c(11, 12), mar = c(0, 0, 5, 0))
for (x in c("Trp63", "Nppc")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:10) plot.new()

for (x in c("Serpinb13")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:11) plot.new()

for (x in c("Mmp10", "Ctgf", "Phlda1", "Cyr61", "Tnfrsf12a")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:7) plot.new()

for (x in c("Tuba1b", "Hist1h4c", "Stmn1", "Pttg1", "Ube2c")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:7) plot.new()

for (x in c("Serpinb3b", "Cldn4", "Krt16", "Plac8", "Plaur")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:7) plot.new()

for (x in c("Jun", "Zfp36", "Junb", "Id1", "Ier2")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:7) plot.new()

for (x in c("Postn", "Islr", "Snca", "Pcdh7", "Adh7", "Alox15")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:6) plot.new()

for (x in c("Mki67", "Top2a", "Pbk", "Gtse1")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:8) plot.new()

for (x in c("Krt7", "Scgb3a2", "Hes1", "Clca2", "Serpinb3a", "Serpinb3b", "Serpinb3c", "Serpinb3d")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:4) plot.new()

for (x in c("Aqp3", "Icam1", "Krt17", "Krt5", "Krt15", "Sfn",
            "Perp", "Fxyd3", "Sdc1", "Gstm2", "F3", "Epas1")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}

for (x in c("Krt5", "Krt8", "Col17a1", "Igfbp3", "Krt14", "Bcam", "Dcn")){
  colors <- tpm[x, clusters %in% 1:6]
  if (sum(colors > 0) == 0){
    colors <- rep("Grey", length(colors))
  }else{
    colors <- colorRampPalette(c("Grey", "Yellow", "Red"))(102)[ceiling((colors-min(colors))/(max(colors)-min(colors))*100)+1]}
  plot(tsne[clusters %in% 1:6, ], pch = 16, cex = 0.5, col = colors, axes = F, xlab = "", ylab = "", main = x, cex.main = 2, font.main = 4)}
for (i in 1:5) plot.new()
dev.off()
#SFigure 5
