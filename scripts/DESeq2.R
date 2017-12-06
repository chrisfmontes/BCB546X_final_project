# Loading required libraries
if (!require("DESeq2")) install.packages("DESeq2") # install DESeq2 if it's not already installed
library("DESeq2")
if (!require("gplots")) install.packages("gplots")
library("gplots")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library("RColorBrewer")

genotypes <- c("Mo17", "B73", "Oh43") #name your different genotypes
treatment <- c("cold", "control", "heat", "Salt" ,"UV") # list your treatments
bioreps <- c(3,3,3,1,1) # write the number of replicates for each treatment, following the order of the "tratment" object
samplenames <- paste(rep(genotypes, each = sum(bioreps)), rep(treatment, bioreps), sep = "_")
counts <- "../E-MTAB-4258-raw-counts.tsv" #here goes the name of your counts table
sampleorder <- read.table("../sample_summary_order.csv", sep = ",", stringsAsFactors = F) # if you want to re-sort the columns, provide a comma-separated list with the intended new order
sampleorder <- unlist(sampleorder[1,])

mircounts <- read.table(file = counts, sep = "\t", header = T) #check whether your count table have a header
row.names(mircounts) <- mircounts[,1]
mircounts <- mircounts[,-c(1,2)]
mircounts <- mircounts[,sampleorder] # skip this if you do not want to re-sort the columns

hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
X11()
barplot(apply(mircounts,2,sum),col=as.factor(samplenames),las=2)
X11()
pairs(log2(mircounts[,c(1:11)]+1),main="Pair-wise sample to sample counts Mo17",pch=".")
X11()
pairs(log2(mircounts[,c(12:22)]+1),main="Pair-wise sample to sample counts B73",pch=".")
X11()
pairs(log2(mircounts[,c(23:33)]+1),main="Pair-wise sample to sample counts Oh43",pch=".")
X11()
pca <- princomp(mircounts)
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="Sample to Sample")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.8)
X11()
heatmap.2(cor(mircounts),trace="none",col=hmcol,main="Sample Correlation")
conds = data.frame(samplenames)
colnames(conds)="treatment"
cds <- DESeqDataSetFromMatrix(countData = mircounts, colData = conds, design = ~ treatment)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- nbinomWaldTest(cds)
X11()
plotDispEsts(cds)
X11()
par(mfrow=c(2,1))
prenorm=apply(mircounts,2,sum)
barplot(prenorm,col=as.factor(samplenames),las=2,names=samplenames)
postnorm=apply(counts(cds,normalized=TRUE),2,sum)
barplot(postnorm,col=as.factor(samplenames),las=2,names=samplenames)
X11()
pca <- princomp(counts(cds,normalized=T))
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="Sample to Sample PCA")
text(pca$loadings, as.vector(colnames(mircounts)), pos=3, cex=0.8)

#we can use the negative-binomial test for each pairwise comparison of interest.
B73_cold =  results( cds, contrast=c("treatment", "B73_cold", "B73_control"))
Mo17_cold <- results( cds, contrast=c("treatment", "Mo17_cold", "Mo17_control"))
Oh43_cold <- results( cds, contrast=c("treatment", "Oh43_cold", "Oh43_control"))

B73_heat =  results( cds, contrast=c("treatment", "B73_heat", "B73_control"))
Mo17_heat <- results( cds, contrast=c("treatment", "Mo17_heat", "Mo17_control"))
Oh43_heat <- results( cds, contrast=c("treatment", "Oh43_heat", "Oh43_control"))

B73_Salt =  results( cds, contrast=c("treatment", "B73_Salt", "B73_control"))
Mo17_Salt <- results( cds, contrast=c("treatment", "Mo17_Salt", "Mo17_control"))
Oh43_Salt <- results( cds, contrast=c("treatment", "Oh43_Salt", "Oh43_control"))

B73_UV =  results( cds, contrast=c("treatment", "B73_UV", "B73_control"))
Mo17_UV <- results( cds, contrast=c("treatment", "Mo17_UV", "Mo17_control"))
Oh43_UV <- results( cds, contrast=c("treatment", "Oh43_UV", "Oh43_control"))

# Sort each result on Adjusted P-Value
B73_cold <- B73_cold[order(B73_cold$padj),]
Mo17_cold <- Mo17_cold[order(Mo17_cold$padj),]
Oh43_cold <- Oh43_cold[order(Oh43_cold$padj),]

B73_heat <- B73_heat[order(B73_heat$padj),]
Mo17_heat <- Mo17_heat[order(Mo17_heat$padj),]
Oh43_heat <- Oh43_heat[order(Oh43_heat$padj),]

B73_Salt <- B73_Salt[order(B73_Salt$padj),]
Mo17_Salt <- Mo17_Salt[order(Mo17_Salt$padj),]
Oh43_Salt <- Oh43_Salt[order(Oh43_Salt$padj),]

B73_UV <- B73_UV[order(B73_UV$padj),]
Mo17_UV <- Mo17_UV[order(Mo17_UV$padj),]
Oh43_UV <- Oh43_UV[order(Oh43_UV$padj),]
# Choose significant miRs for each contrast by log fold change and adj. P-value
B73_cold_sig = rownames(B73_cold[(abs(B73_cold$log2FoldChange) > 2) & (B73_cold$padj < 0.1) & !is.na(B73_cold$padj),])
siglist_cold = unique(B73_cold_sig)
Mo17_cold_sig = rownames(Mo17_cold[(abs(Mo17_cold$log2FoldChange) > 2) & (Mo17_cold$padj < 0.1) & !is.na(Mo17_cold$padj),])
siglist2_cold = unique(Mo17_cold_sig)
Oh43_cold_sig = rownames(Oh43_cold[(abs(Oh43_cold$log2FoldChange) > 2) & (Oh43_cold$padj < 0.1) & !is.na(Oh43_cold$padj),])
siglist3_cold = unique(Oh43_cold_sig)

B73_heat_sig = rownames(B73_heat[(abs(B73_heat$log2FoldChange) > 2) & (B73_heat$padj < 0.1) & !is.na(B73_heat$padj),])
siglist_heat = unique(B73_heat_sig)
Mo17_heat_sig = rownames(Mo17_heat[(abs(Mo17_heat$log2FoldChange) > 2) & (Mo17_heat$padj < 0.1) & !is.na(Mo17_heat$padj),])
siglist2_heat = unique(Mo17_heat_sig)
Oh43_heat_sig = rownames(Oh43_heat[(abs(Oh43_heat$log2FoldChange) > 2) & (Oh43_heat$padj < 0.1) & !is.na(Oh43_heat$padj),])
siglist3_heat = unique(Oh43_heat_sig)

B73_Salt_sig = rownames(B73_Salt[(abs(B73_Salt$log2FoldChange) > 2) & (B73_Salt$padj < 0.1) & !is.na(B73_Salt$padj),])
siglist_Salt = unique(B73_Salt_sig)
Mo17_Salt_sig = rownames(Mo17_Salt[(abs(Mo17_Salt$log2FoldChange) > 2) & (Mo17_Salt$padj < 0.1) & !is.na(Mo17_Salt$padj),])
siglist2_Salt = unique(Mo17_Salt_sig)
Oh43_Salt_sig = rownames(Oh43_Salt[(abs(Oh43_Salt$log2FoldChange) > 2) & (Oh43_Salt$padj < 0.1) & !is.na(Oh43_Salt$padj),])
siglist3_Salt = unique(Oh43_Salt_sig)

B73_UV_sig = rownames(B73_UV[(abs(B73_UV$log2FoldChange) > 2) & (B73_UV$padj < 0.1) & !is.na(B73_UV$padj),])
siglist_UV = unique(B73_UV_sig)
Mo17_UV_sig = rownames(Mo17_UV[(abs(Mo17_UV$log2FoldChange) > 2) & (Mo17_UV$padj < 0.1) & !is.na(Mo17_UV$padj),])
siglist2_UV = unique(Mo17_UV_sig)
Oh43_UV_sig = rownames(Oh43_UV[(abs(Oh43_UV$log2FoldChange) > 2) & (Oh43_UV$padj < 0.1) & !is.na(Oh43_UV$padj),])
siglist3_UV = unique(Oh43_UV_sig)

siglist_tot <- unique(c(siglist_cold, siglist_heat, siglist_Salt, siglist_UV))

#Lets make some volcano plots of the comparison
X11()
plot(B73_cold$log2FoldChange,-log(B73_cold$padj,10),main="Volcano Plot B73 cold treated vs control")
text(B73_cold[1:20,]$log2FoldChange,-log(B73_cold[1:20,]$padj,10),labels=rownames(B73_cold[1:20,]),cex=0.7,pos=1)
legend("topleft","control",cex=0.5)
legend("topright","induced",cex=0.5)

# Generate Nice heatmap colours
hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
# Heatmap of significant hits
X11()
heatmap.2(log2(counts(cds[siglist_tot,c(1:3,7:14,18:25,29:33)],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames[c(1:3,7:14,18:25,29:33)],margin=c(7,10),main = "DE genes per treatment")
png(filename = "../heat.png")
heatmap.2(log2(counts(cds[siglist_tot,c(1:3,7:14,18:25,29:33)],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames[c(1:3,7:14,18:25,29:33)],margin=c(7,10),main = "DE genes per treatment")
dev.off()
x11()
heatmap.2(log2(counts(cds[siglist2,],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames,margin=c(7,10),main = "log2FoldChange > 1.0")
#finally we write a file with all the genes and a pdf with the graphics
write.table(B73_cold,file = "DESeq2.table.txt", sep = "\t", quote = FALSE)
pdf(file = "DESeq2.graphics.pdf")
plot(B73_cold$log2FoldChange,-log(B73_cold$padj,10),main="Volcano Plot control vs induced")
text(B73_cold[1:20,]$log2FoldChange,-log(B73_cold[1:20,]$padj,10),labels=rownames(B73_cold[1:20,]),cex=0.7,pos=1)
legend("topleft","control",cex=0.5)
legend("topright","induced",cex=0.5)
heatmap.2(log2(counts(cds[siglist,],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames,margin=c(7,10),main = "log2FoldChange > 1.5")
heatmap.2(log2(counts(cds[siglist2,],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames,margin=c(7,10),main = "log2FoldChange > 1.0", cexRow = 0.3)
dev.off()