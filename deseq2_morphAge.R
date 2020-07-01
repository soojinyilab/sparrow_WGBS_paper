library(DESeq2)

### Read data
cond <- "female_Hyp"
cts <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id", check.names = F)) # raw count
coldata <- read.csv(paste(cond, ".cond", sep = ""), sep="\t", row.names=1, header=F) # sample condition file
colnames(coldata) <- c("morph", "age", "type")

rownames(coldata) <- sub("fb", "", rownames(coldata))
print(all(rownames(coldata) %in% colnames(cts)))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

### Differential expression between two morphs
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ age + morph) # corrects for age
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, fitType='local')
counts <- counts(dds, normalized=TRUE)

res <- results(dds)
baseMeanSep <- sapply(levels(dds$morph), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$morph == lvl]))
res <- cbind(as.data.frame(res), baseMeanSep)
res$Gene <- rownames(res)
res <- res[,c(9,1,7,8,2,3,4,5,6)]
res <- cbind(res, counts)

resOrdered <- res[order(res$padj), ]
write.table(resOrdered, file=paste(cond, ".morph.DESeq2", sep = ""), quote=F, sep="\t", row.names=F)
write.table(cbind(dds$morph, dds$sizeFactor), file=paste(cond, "_sizeFactor.txt", sep=""), quote=F, sep="\t", col.names=F)

## Differential expression between two age groups
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ morph + age) # corrects for morph
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, fitType='local')
counts <- counts(dds, normalized=TRUE)

res <- results(dds)
baseMeanSep <- sapply(levels(dds$age), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$age == lvl]))
res <- cbind(as.data.frame(res), baseMeanSep)
res$Gene <- rownames(res)
res <- res[,c(9,1,7,8,2,3,4,5,6)]
res <- cbind(res, counts)

resOrdered <- res[order(res$padj), ]
write.table(resOrdered, file=paste(cond, ".age.DESeq2", sep = ""), quote=F, sep="\t", row.names=F)


