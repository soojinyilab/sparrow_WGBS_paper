library(DSS)

### Read formatted bismark coverage file into a DSS object 
samples <- c('N1-C10-01_ZAL2','N1-C10-08_ZAL2','N3-C10-19_ZAL2',
             '23614-WSF_ZAL2','WSF-13-02_ZAL2','WSF-13-03_ZAL2','WSF-13-08_ZAL2',
             'N1-C10-01_ZAL2m','N1-C10-08_ZAL2m','N3-C10-19_ZAL2m',
             '23614-WSF_ZAL2m','WSF-13-02_ZAL2m','WSF-13-03_ZAL2m','WSF-13-08_ZAL2m')

data <- list()
for (i in 1:length(samples)) {
  print(samples[i])
  data[[i]] = read.table(paste(samples[i], ".mCount.txt", sep=""), col.names = c("chr","pos","N","X")) # see DSS documentation for the format requirement
}

BSdat <- makeBSseqData(data, samples) #  CpGs
BSdat <- sort(BSdat)

### Make a design table
design <- data.frame(chrom = c(rep("WS_ZAL2",7),
                               rep("WS_ZAL2m",7)),
                     age = c(rep("Chick",3), rep("Adult",4),
                             rep("Chick",3), rep("Adult",4)))
design$chrom <- factor(design$chrom, levels = c("WS_ZAL2","WS_ZAL2m"))
design$age <- factor(design$age, levels = c("Adult","Chick"))
row.names(design) <- samples
pData(BSdat) <- design

### Filter CpGs
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(BSdat, type="Cov") == 0) == 0)

length(loci.idx)
BSdat.filtered <- BSdat[loci.idx,]

##### Identify DMLs #####
design$chrom <- factor(design$chrom, levels = c("WS_ZAL2","WS_ZAL2m"))
DMLfit = DMLfit.multiFactor(BSdat.filtered, design=design, formula=~chrom+age) # WS_ZAL2 as the reference level
coefs <- colnames(DMLfit$X)
statL <- list() # list of stats
pvalsL <- list() # list of pvalues
fdrsL <- list() # list of qvalues
bonfL <- list()
for (coef in coefs[2:length(coefs)]) {
  DMLtest = DMLtest.multiFactor(DMLfit, coef=coef) # DMLfit$X can check the term to test
  statL[[paste0('stat_',coef)]] = DMLtest$stat
  pvalsL[[paste0('pval_',coef)]] = DMLtest$pvals
  fdrsL[[paste0('fdr_',coef)]] = DMLtest$fdrs
  bonfL[[paste0('bonf_',coef)]] = p.adjust(DMLtest$pvals, method = "bonferroni")
}

### Output all CpGs
DMLtest = DMLtest.multiFactor(DMLfit, coef="chromWS_ZAL2m")

methy.df <- data.frame(DMLtest)
methy.df$start <- methy.df$pos
methy.df$end <- methy.df$pos
methy.df <- subset(methy.df, select=c("chr", "start", "end"))

WS_ZAL2 <- rownames(subset(design, chrom == "WS_ZAL2"))
WS_ZAL2m <- rownames(subset(design, chrom == "WS_ZAL2m"))
gps <- c("WS_ZAL2","WS_ZAL2m")

for (gp in gps) {
  print(gp)
  print(get(gp))
  methy.df[,gp] <- apply(getMeth(BSdat.filtered, region=methy.df, type = "raw", what = "perRegion")[,get(gp)], FUN=mean, MARGIN=1,na.rm=TRUE)
}

methy.df <- cbind(methy.df,
                  diff_chromWS_ZAL2m = methy.df$WS_ZAL2m - methy.df$WS_ZAL2)
methyDat <- getMeth(BSdat.filtered, region=methy.df, type = "raw", what = "perRegion")[,samples]
methy.df <- cbind(methy.df, as.data.frame(statL), as.data.frame(pvalsL), as.data.frame(fdrsL), as.data.frame(bonfL), methyDat)

write.table(methy.df, "methyDSS.txt", row.names = F, quote = F, sep = "\t")
saveRDS(methy.df, "methyDF.rds")
save.image("dss_allele.RData")
