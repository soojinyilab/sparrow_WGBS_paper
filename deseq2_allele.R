library(DESeq2)

### Read data
cond <- "female_Hyp"
cts <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id", check.names = F)) # raw count
coldata <- read.csv(paste(cond, ".cond", sep = ""), sep="\t", row.names=1, header=F) # sample condition file
colnames(coldata) <- c("ind","allele","age","type")  
coldata$allele <- factor(coldata$allele, levels = c("ZAL2","ZAL2m"))

rownames(coldata) <- sub("fb", "", rownames(coldata))
print(all(rownames(coldata) %in% colnames(cts)))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

### Differential expression between two alleles
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ age + allele) # corrects for age
sizeFactors.df <- read.table(paste(cond, "_sizeFactor.txt", sep = "")) # size factors should come from deseq2_morphAge.R
sizeFactors(dds) <- rep(subset(sizeFactors.df, V2 == 2)$V3, each = 2)
counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds, fitType='local')
print(as.data.frame(dds@colData))

res <- results(dds)
baseMeanSep <- sapply(levels(dds$allele), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$allele == lvl]))
res <- cbind(as.data.frame(res), baseMeanSep)
res$Gene <- rownames(res)
res <- res[,c(9,1,7,8,2,3,4,5,6)]
res <- cbind(res, counts)

resOrdered <- res[order(res$padj), ]
write.table(resOrdered, file=paste(cond, ".allele.DESeq2", sep = ""), quote=F, sep="\t", row.names=F)

