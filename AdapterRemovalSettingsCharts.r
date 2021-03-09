#!/usr/bin/env Rscript
#Code by Lindsey Fenderson - February 2021
args = commandArgs(trailingOnly = TRUE)
#
args <- commandArgs()
print(args)
FILENAME <- args[6]
print (FILENAME)

TrimmingStats <- read.table(FILENAME, header = TRUE)

options(scipen=10)
pdf(file=paste(FILENAME, "All.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$All, type="l", col="navy", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="All Reads", main="Number and Size of All Reads \nBefore Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length))) #cex.axis=0.35)
axis(side=2, 
     at = pretty(range(TrimmingStats$All)),las=2)
mtext(text = "Read Length",
      side = 1,#side 1 = bottom
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "Mate1.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$Mate1, type="l", col="darkgreen", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Mate1", main="Number and Size of Retained Mate1 Reads \nFollowing Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$Mate1)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "Mate2.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$Mate2, type="l", col="aquamarine1", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Mate2", main="Number and Size of Retained Mate2 Reads \nFollowing Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$Mate2)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "Singletons.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$Singleton, type="l", col="purple4", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Singletons", main="Number and Size of Retained Singleton Reads \nFollowing Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$Singleton)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "Collapsed.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$Collapsed, type="l", col="violetred", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Collapsed", main="Number and Size of Retained Collapsed Reads \nFollowing Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$Collapsed)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "CollapsedTruncated.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$CollapsedTruncated, type="l", col="cadetblue1", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Collapsed & Truncated", main="Number and Size of Retained Collapsed & Truncated \nReads Following Quality Trimming by Adapter Removal", cex.main=0.9)
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$CollapsedTruncated)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

pdf(file=paste(FILENAME, "Discarded.pdf", sep="-"), width = 6, height = 4, family = "NewCenturySchoolbook")
par(mar = c(4.5, 6.5, 2.5, 2.5), mgp = c(5, 1, 0))
plot(TrimmingStats$Length,TrimmingStats$Discarded, type="l", col="black", lwd=3, xaxt="none", yaxt="none", xlab="", ylab="Mate1", main="Number and Size of Discarded Reads \nFollowing Quality Trimming by Adapter Removal")
axis(side=1,
     at = pretty(range(TrimmingStats$Length)))
axis(side=2,
     at = pretty(range(TrimmingStats$Discarded)),las=2)
mtext(text = "Read Length",
      side = 1,
      line = 3)
dev.off()

q()
n

