## Begin of the DESeq2 analysis script ##

## Loading required libraries ##

# Load DESeq2, gplots and RColorBrewer libraries if present, install them if not present.
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("gplots")) install.packages("gplots")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

## Setting up the data structure for DESeq ##

# First steps: defining variales and pathways
genotypes <- c("Mo17", "B73", "Oh43") #name your different genotypes
treatment <- c("cold", "control", "heat", "Salt" ,"UV") # list your treatments
bioreps <- c(3,3,3,1,1) # write the number of replicates for each treatment, following the order of the "tratment" object
samplenames <- paste(rep(genotypes, each = sum(bioreps)), rep(treatment, bioreps), sep = "_")
counts <- "../tables_and_data/counts.txt" #here goes the path and name of your counts table
sampleorder <- read.table("../tables_and_data/sample_summary_order.csv", sep = ",", stringsAsFactors = F) # if you want to re-sort the columns, provide a comma-separated list with the intended new order
sampleorder <- unlist(sampleorder[1,])

# Second steps: setting up de counts table
mycounts <- read.table(file = counts, sep = " ", header = T) #check whether your count table have a header and the field delimiter
row.names(mycounts) <- mycounts[,1]
mycounts <- mycounts[,-1]
mycounts <- mycounts[,sampleorder] # skip this if you do not want to re-sort the columns

# Third step (optional): creating some diagnostic plots
hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(100) # set some nice colors
# Bar plot of counts for each sample
X11()
par(mar=c(6.5,4.5,4.1,2.1))
barplot(apply(mycounts,2,sum),col=as.factor(samplenames),las=2,names=samplenames)
# Pairwise scatterplot for all samples in each gentype
X11()
pairs(log2(mycounts[,c(1:11)]+1),main="Pair-wise sample to sample counts Mo17",pch=".")
X11()
pairs(log2(mycounts[,c(12:22)]+1),main="Pair-wise sample to sample counts B73",pch=".")
X11()
pairs(log2(mycounts[,c(23:33)]+1),main="Pair-wise sample to sample counts Oh43",pch=".")
# Principal component analysis plot for all samples
X11()
pca <- princomp(mycounts)
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="Sample to Sample")
text(pca$loadings, as.vector(colnames(mycounts)), pos=3, cex=0.8)
# We can now export all these plot as png image files
# Exporting the pairwise scatterplots
png(filename = "../scatter_Mo17.png", width = 1000, height = 1000)
pairs(log2(mycounts[,c(1:11)]+1),main="Pair-wise sample to sample counts Mo17",pch=".", cex.axis = 2.0, cex.labels = 1.5, cex.main = 2.0)
dev.off()
png(filename = "../scatter_B73.png", width = 1000, height = 1000)
pairs(log2(mycounts[,c(12:22)]+1),main="Pair-wise sample to sample counts B73",pch=".", cex.axis = 2.0, cex.labels = 1.5, cex.main = 2.0)
dev.off()
png(filename = "../scatter_Oh43.png", width = 1000, height = 1000)
pairs(log2(mycounts[,c(23:33)]+1),main="Pair-wise sample to sample counts Oh43",pch=".", cex.axis = 2.0, cex.labels = 1.5, cex.main = 2.0)
dev.off()
# Exporting the PCA plot
png(filename = "../PCA_prenorm.png", width = 1000, height = 1000)
pca <- princomp(mycounts)
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="Pre-normalization PCA", cex.main = 2)
text(pca$loadings, as.vector(colnames(mycounts)), pos=3, cex=0.8)
dev.off()

## DESeq2 analysis
# For this analysis we create a conds object that states the treatments to be used for DE profile
conds = data.frame(samplenames)
colnames(conds)="treatment"
# We run the DESeq2 normalization and create a DESeq type object called cds (it will contain all the DE analysis data)
cds <- DESeqDataSetFromMatrix(countData = mycounts, colData = conds, design = ~ treatment)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- nbinomWaldTest(cds)

# Now we can run some diagnostic plots again
# Bar plot of counts for each sample, comparing before and after normalization
X11()
par(mfrow=c(2,1), mar=c(6.5,4.5,4.1,2.1))
prenorm=apply(mycounts,2,sum)
barplot(prenorm,col=as.factor(samplenames),las=2,names=samplenames)
postnorm=apply(counts(cds,normalized=TRUE),2,sum)
barplot(postnorm,col=as.factor(samplenames),las=2,names=samplenames)
# Principal component analysis after normalization
X11()
pca <- princomp(counts(cds,normalized=T))
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="PCA post-norm")
text(pca$loadings, as.vector(colnames(mycounts)), pos=3, cex=0.8)
# We now export these plots as png image files
# Exporting the bar plot, comparing before and after normalization
png(filename = "../samp_barplot_comp.png", width = 1000, height = 1000)
par(mfrow=c(2,1), mar=c(6.5,4.5,4.1,2.1))
prenorm=apply(mycounts,2,sum)
barplot(prenorm,col=as.factor(samplenames),las=2,names=samplenames)
postnorm=apply(counts(cds,normalized=TRUE),2,sum)
barplot(postnorm,col=as.factor(samplenames),las=2,names=samplenames)
dev.off()
# Exporting the PCA plot
png(filename = "../PCA_postnorm.png", width = 1000, height = 1000)
pca <- princomp(counts(cds,normalized=T))
plot(pca$loadings, col=as.factor(samplenames),  pch=19, cex=2, main="Post-normalization PCA", cex.main = 2)
text(pca$loadings, as.vector(colnames(mycounts)), pos=3, cex=0.8)
dev.off()

# We can use the negative-binomial test for each pairwise comparison of interest.
# Cold treated samples vs control
B73_cold =  results( cds, contrast=c("treatment", "B73_cold", "B73_control"))
Mo17_cold <- results( cds, contrast=c("treatment", "Mo17_cold", "Mo17_control"))
Oh43_cold <- results( cds, contrast=c("treatment", "Oh43_cold", "Oh43_control"))

# Heat-treated samples vs control
B73_heat =  results( cds, contrast=c("treatment", "B73_heat", "B73_control"))
Mo17_heat <- results( cds, contrast=c("treatment", "Mo17_heat", "Mo17_control"))
Oh43_heat <- results( cds, contrast=c("treatment", "Oh43_heat", "Oh43_control"))

# Salt-treated samples vs control
B73_Salt =  results( cds, contrast=c("treatment", "B73_Salt", "B73_control"))
Mo17_Salt <- results( cds, contrast=c("treatment", "Mo17_Salt", "Mo17_control"))
Oh43_Salt <- results( cds, contrast=c("treatment", "Oh43_Salt", "Oh43_control"))

# UV-treated samples vs control
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

# Choose significantly DE genes for each contrast by log fold change and adj. P-value
# log2 fold over 2 and padj less 0.1

# Cold treatment DE genes
B73_cold_sig = rownames(B73_cold[(abs(B73_cold$log2FoldChange) > 2) & (B73_cold$padj < 0.1) & !is.na(B73_cold$padj),])
siglist_cold = unique(B73_cold_sig)
Mo17_cold_sig = rownames(Mo17_cold[(abs(Mo17_cold$log2FoldChange) > 2) & (Mo17_cold$padj < 0.1) & !is.na(Mo17_cold$padj),])
siglist2_cold = unique(Mo17_cold_sig)
Oh43_cold_sig = rownames(Oh43_cold[(abs(Oh43_cold$log2FoldChange) > 2) & (Oh43_cold$padj < 0.1) & !is.na(Oh43_cold$padj),])
siglist3_cold = unique(Oh43_cold_sig)

# Heat treatment DE genes
B73_heat_sig = rownames(B73_heat[(abs(B73_heat$log2FoldChange) > 2) & (B73_heat$padj < 0.1) & !is.na(B73_heat$padj),])
siglist_heat = unique(B73_heat_sig)
Mo17_heat_sig = rownames(Mo17_heat[(abs(Mo17_heat$log2FoldChange) > 2) & (Mo17_heat$padj < 0.1) & !is.na(Mo17_heat$padj),])
siglist2_heat = unique(Mo17_heat_sig)
Oh43_heat_sig = rownames(Oh43_heat[(abs(Oh43_heat$log2FoldChange) > 2) & (Oh43_heat$padj < 0.1) & !is.na(Oh43_heat$padj),])
siglist3_heat = unique(Oh43_heat_sig)

# Salt treatment DE genes
B73_Salt_sig = rownames(B73_Salt[(abs(B73_Salt$log2FoldChange) > 2) & (B73_Salt$padj < 0.1) & !is.na(B73_Salt$padj),])
siglist_Salt = unique(B73_Salt_sig)
Mo17_Salt_sig = rownames(Mo17_Salt[(abs(Mo17_Salt$log2FoldChange) > 2) & (Mo17_Salt$padj < 0.1) & !is.na(Mo17_Salt$padj),])
siglist2_Salt = unique(Mo17_Salt_sig)
Oh43_Salt_sig = rownames(Oh43_Salt[(abs(Oh43_Salt$log2FoldChange) > 2) & (Oh43_Salt$padj < 0.1) & !is.na(Oh43_Salt$padj),])
siglist3_Salt = unique(Oh43_Salt_sig)

# UV treatment DE genes
B73_UV_sig = rownames(B73_UV[(abs(B73_UV$log2FoldChange) > 2) & (B73_UV$padj < 0.1) & !is.na(B73_UV$padj),])
siglist_UV = unique(B73_UV_sig)
Mo17_UV_sig = rownames(Mo17_UV[(abs(Mo17_UV$log2FoldChange) > 2) & (Mo17_UV$padj < 0.1) & !is.na(Mo17_UV$padj),])
siglist2_UV = unique(Mo17_UV_sig)
Oh43_UV_sig = rownames(Oh43_UV[(abs(Oh43_UV$log2FoldChange) > 2) & (Oh43_UV$padj < 0.1) & !is.na(Oh43_UV$padj),])
siglist3_UV = unique(Oh43_UV_sig)

siglist_tot <- unique(c(siglist_cold, siglist_heat, siglist_Salt, siglist_UV))

# Lets make some cool plots of the comparison
# Volcano plots of cold treated plants
X11()
par(mfrow=c(1,3))
plot(B73_cold$log2FoldChange,-log(B73_cold$padj,10),main="B73 cold treated vs control", ylab = "-log(padj)", xlab = "log2fold (cold/control)")
text(B73_cold[1:20,]$log2FoldChange,-log(B73_cold[1:20,]$padj,10), labels=rownames(B73_cold[1:20,]),cex=0.7,pos=1)
plot(Mo17_cold$log2FoldChange,-log(Mo17_cold$padj,10),main="Mo17 cold treated vs control", ylab = "-log(padj)", xlab = "log2fold (cold/control)")
text(Mo17_cold[1:20,]$log2FoldChange,-log(Mo17_cold[1:20,]$padj,10), labels=rownames(Mo17_cold[1:20,]),cex=0.7,pos=1)
plot(Oh43_cold$log2FoldChange,-log(Oh43_cold$padj,10),main="Oh43 cold treated vs control", ylab = "-log(padj)", xlab = "log2fold (cold/control)")
text(Oh43_cold[1:20,]$log2FoldChange,-log(Oh43_cold[1:20,]$padj,10), labels=rownames(Oh43_cold[1:20,]),cex=0.7,pos=1)
# Generate Nice heatmap colours (if not loaded already)
hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
# Heatmap of significant hits
X11()
heatmap.2(log2(counts(cds[siglist_tot,c(1:3,7:14,18:25,29:33)],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames[c(1:3,7:14,18:25,29:33)],margin=c(7,10),main = "DE genes per treatment")
# Export heatmap to a png image file
png(filename = "../heat.png", width = 1000, height = 1000)
heatmap.2(log2(counts(cds[siglist_tot,c(1:3,7:14,18:25,29:33)],normalized=TRUE)+1),col=hmcol,trace="none",labCol=samplenames[c(1:3,7:14,18:25,29:33)],margins=c(7,10),main = "DE genes per treatment", cexCol = 1.4)
dev.off()

#finally we write a file with all the genes for each condition and treatment
write.table(B73_cold,file = "../B73_cold_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Mo17_cold,file = "../Mo17_cold_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Oh43_cold,file = "../Oh43_cold_DE.genes_table.txt", sep = "\t", quote = FALSE)

write.table(B73_heat,file = "../B73_heat_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Mo17_heat,file = "../Mo17_heat_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Oh43_heat,file = "../Oh43_heat_DE.genes_table.txt", sep = "\t", quote = FALSE)

write.table(B73_Salt,file = "../B73_Salt_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Mo17_Salt,file = "../Mo17_Salt_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Oh43_Salt,file = "../Oh43_Salt_DE.genes_table.txt", sep = "\t", quote = FALSE)

write.table(B73_UV,file = "../B73_UV_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Mo17_UV,file = "../Mo17_UV_DE.genes_table.txt", sep = "\t", quote = FALSE)
write.table(Oh43_UV,file = "../Oh43_UV_DE.genes_table.txt", sep = "\t", quote = FALSE)