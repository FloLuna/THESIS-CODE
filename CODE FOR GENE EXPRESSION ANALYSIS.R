# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)

# load series and platform data from GEO

gset <- getGEO("GSE95404", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("CSF1","CSF2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=1000000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC",
                          "Platform_SPOTID",
                          "Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE95404", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE95404", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE95404")
ex_df <- as.data.frame(ex)
ex_df <- rownames_to_column(ex_df, var = 'ID')
ex_df$ID <- as.integer(ex_df$ID)
combined1 <- inner_join(tT, ex_df, by = 'ID')
combined1 <- inner_join(combined1, tT2, by = 'ID')





library(ggplot2)



#-----------------------------------------------
#.......CHOOSE YOUR OWN GENES HERE
#______________________________________________


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Remove leading/trailing whitespaces from sample names
group_data$Sample <- trimws(group_data$Sample)
colnames(combined1) <- trimws(colnames(combined1))

# Reshape the data into a long format
long_data <- combined1 %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Remove leading/trailing whitespaces from sample names
long_data$Sample <- trimws(long_data$Sample)

# Match samples with groups using the match function
long_data$Group <- group_data$Group[match(long_data$Sample, group_data$Sample)]

selected_genes <- c('Il6', "Nlrp3", "Cd86", "Ifitm1", 'Itgax', 'Apobec1', 'Tlr8', 'Oas2', 'Tmem173', 'Gbp5', 'Tlr7',
                    'Samhd1', 'Oasl2', 'Ifih1', 'Bnip3l', 'Atg7', 'Plscr1', 'Oas3',
                    'Itch', 'Ptprc', 'Dhx58', 'Cd40', 'Ddx41', 'Rsad2', 'Irf5', 'Ctsl')

subset_data <- long_data %>%
  filter(Gene.symbol.x %in% selected_genes)

# Calculate the absolute difference between average logFC values of two groups
ordered_data <- subset_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)

# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))


#-------------------------------------------
#........CYTOKINE DATA HERE
#____________________________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"

cytokine_data <- combined1 %>%
  filter(grepl(cytokine_ID, GO.Function.ID))
  
cytokine_data <- cytokine_data %>%
  filter(!grepl(chemokine_ID, GO.Function.ID))

library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
cytokine_long_data <- cytokine_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
cytokine_long_data$Group <- group_data$Group[match(cytokine_long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- cytokine_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[150:186]


selected_genes <- c(top_25_genes)

selected_genes

subset_data <- cytokine_long_data %>%
  filter(Gene.symbol.x %in% selected_genes)

# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#--------------------------------------------------------
#........ CHEMOKINE DATA HERE
#________________________________________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"

chemokine_data <- combined1 %>%
  filter(grepl(chemokine_ID, GO.Function.ID))


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
chemokine_long_data <- chemokine_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
chemokine_long_data$Group <- group_data$Group[match(chemokine_long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- chemokine_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[30:36]


selected_genes <- c(top_25_genes)

selected_genes

subset_data <- chemokine_long_data %>%
  filter(Gene.symbol.x %in% selected_genes)


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#---------------------------------------------
#............CYTOKINE RECEPTOR DATA HERE
#____________________________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"
cytokinereceptor_ID <- '4896'
chemokinereceptor_ID <- '4950'

cytokinereceptor_data <- combined1 %>%
  filter(grepl(cytokinereceptor_ID, GO.Function.ID))


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
cytokinereceptor_long_data <- cytokinereceptor_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
cytokinereceptor_long_data$Group <- group_data$Group[match(cytokinereceptor_long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- cytokinereceptor_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[20:42]

selected_genes <- c(top_25_genes)

selected_genes

subset_data <- cytokinereceptor_long_data %>%
  filter(Gene.symbol.x %in% selected_genes)


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#--------------------------------------------------
#..................CHEMOKINE RECEPTOR DATA HERE
#__________________________________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"
cytokinereceptor_ID <- '4896'
chemokinereceptor_ID <- '4950'

chemokinereceptor_data <- combined1 %>%
  filter(grepl(chemokinereceptor_ID, GO.Function.ID))


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
chemokinereceptor_long_data <- chemokinereceptor_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
chemokinereceptor_long_data$Group <- group_data$Group[match(chemokinereceptor_long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- chemokinereceptor_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[14:19]


selected_genes <- c(top_25_genes)

selected_genes

subset_data <- chemokinereceptor_long_data %>%
  filter(Gene.symbol.x %in% selected_genes)


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#------------------------------------------------
#............. INDIVIDUAL SAMPLES HERE
#_______________________________________________
library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Remove leading/trailing whitespaces from sample names
group_data$Sample <- trimws(group_data$Sample)
colnames(combined1) <- trimws(colnames(combined1))

# Reshape the data into a long format
long_data <- combined1 %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Remove leading/trailing whitespaces from sample names
long_data$Sample <- trimws(long_data$Sample)

# Match samples with groups using the match function
long_data$Group <- group_data$Group[match(long_data$Sample, group_data$Sample)]

selected_genes <- c('Il6', "Nlrp3", "Cd86", "Ifitm1", 'Itgax', 'Apobec1', 'Tlr8', 'Oas2', 'Tmem173', 'Gbp5', 'Tlr7',
                    'Samhd1', 'Oasl2', 'Ifih1', 'Bnip3l', 'Atg7', 'Plscr1', 'Oas3',
                    'Itch', 'Ptprc', 'Dhx58', 'Cd40', 'Ddx41', 'Rsad2', 'Irf5')

subset_data <- long_data %>%
  filter(Gene.symbol.x %in% selected_genes)

# Calculate the absolute difference between average logFC values of two groups
ordered_data <- subset_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)

# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with samples
ggplot(data = subset_data, aes(x = Sample, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(4, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#---------------------------------------------
#..............TOP GENES FROM ALL
#_____________________________________________

library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Remove leading/trailing whitespaces from sample names
group_data$Sample <- trimws(group_data$Sample)
colnames(combined1) <- trimws(colnames(combined1))

# Reshape the data into a long format
long_data <- combined1 %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
long_data$Group <- group_data$Group[match(long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[22075:22104]
top_25_genes

selected_genes <- c(top_25_genes)

selected_genes

subset_data <- long_data %>%
  filter(Gene.symbol.x %in% selected_genes)


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 7,
                       limits = c(3, 11),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#IF YOU WANNA CHECK THE ROLE RUN THIS

combined1$GO.Process[combined1$Gene.symbol.x == "Ddit4"]

#--------------------------------------------
#...........................VIRAL DATA HERE
#______________________________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"
cytokinereceptor_ID <- '4896'
chemokinereceptor_ID <- '4950'

viral_data <- combined1 %>%
  filter(grepl('GO:0051607', GO.Process.ID))


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
viral_long_data <- viral_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
viral_long_data$Group <- group_data$Group[match(viral_long_data$Sample, group_data$Sample)]


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- viral_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[110:139]


selected_genes <- c(top_25_genes)

selected_genes

subset_data <- viral_long_data %>%
  filter(Gene.symbol.x %in% selected_genes)


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                      breaks = c(seq(0, 10, 0.1)),
                      labels = labels,
                      low = 'blue', mid = 'white', high = 'red',
                      midpoint = 6.5,
                      limits = c(3, 10),
                      oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))

#----------------------------------
#...............MISC
#___________________________

library(dplyr)

cytokine_ID <- "5125"
chemokine_ID <- "8009"
cytokinereceptor_ID <- '4896'
chemokinereceptor_ID <- '4950'

viral_data <- combined1 %>%
  filter(grepl('virus', GO.Process))


library(tidyverse)

# Define groups for samples
group_data <- data.frame(Sample = colnames(combined1)[10:15],
                         Group = c("CSF1", "CSF1", "CSF1", "CSF2", "CSF2", "CSF2"))

# Reshape the data into a long format
viral_long_data <- viral_data %>%
  select(Gene.symbol.x, starts_with("GSM")) %>%
  pivot_longer(cols = starts_with("GSM"), names_to = "Sample", values_to = "logFC")

# Match samples with groups using the match function
viral_long_data$Group <- group_data$Group[match(viral_long_data$Sample, group_data$Sample)]

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# Calculate the absolute difference between average logFC values of two groups
ordered_data <- viral_long_data %>%
  group_by(Gene.symbol.x, Group) %>%
  summarize(avg_logFC = mean(logFC)) %>%
  pivot_wider(names_from = Group, values_from = avg_logFC) %>%
  mutate(diff_logFC = abs(CSF1 - CSF2))

# Order the genes based on the absolute difference
ordered_genes <- ordered_data %>%
  arrange(diff_logFC) %>%
  pull(Gene.symbol.x)
ordered_genes

top_25_genes <- ordered_genes[300:331]



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>

selected_genes <- c(top_25_genes)

selected_genes

subset_data <- viral_long_data %>%
  filter(Gene.symbol.x %in% selected_genes) %>%
  filter(Gene.symbol.x != 'Hist1h2al///Hist1h2ao///Hist2h2ab///Hist2h2aa2///Hist1h2ai///Hist2h2ac///Hist1h2af///Hist1h2ab///Hist1h2ap///Hist1h2an///Hist1h2ak///Hist1h2ah///Hist1h2ag///Hist1h2ae///Hist1h2ad///Hist1h2ac///Hist1h2aa///Hist3h2a///H2afj///H2afx///Hist2h2aa1///Hist1h3d
  ') %>%
  filter(Gene.symbol.x != 'Hist1h2al///Hist1h2ao///Hist2h2ab///Hist2h2aa2///Hist1h2ai///Hist2h2ac///Hist1h2af///Hist1h2ab///Hist1h2ap///Hist1h2an///Hist1h2ak///Hist1h2ah///Hist1h2ag///Hist1h2ae///Hist1h2ad///Hist1h2ac///Hist1h2aa///Hist3h2a///H2afj///H2afx///Hist2h2aa1///Hist1h3d')


# Define the order of genes
gene_order <- c(ordered_genes)

subset_data$Gene.symbol.x <- factor(subset_data$Gene.symbol.x, levels = gene_order)

labels <- ifelse(floor(seq(0, 10, 0.1)) == seq(0, 10), seq(0, 10), "")

# Create the heatmap with grouped samples
ggplot(data = subset_data, aes(x = Group, y = Gene.symbol.x, fill = logFC)) +
  geom_tile(color = 'black', size = 1) +
  scale_fill_gradient2(name = "Log2 Fold Change",
                       breaks = c(seq(0, 10, 0.1)),
                       labels = labels,
                       low = 'blue', mid = 'white', high = 'red',
                       midpoint = 6.5,
                       limits = c(3, 10),
                       oob = scales::squish) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(legend.key.height = unit(2, "cm"))


