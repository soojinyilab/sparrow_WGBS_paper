library(DSS)

##### Read formatted bismark coverage file into a DSS object #####
samples <- c('N1-C10-01','N1-C10-08','N3-C10-19',
             '23614-WSF','WSF-13-02','WSF-13-03','WSF-13-08',
             'N1-C10-10','N1-C10-19',
             'TSF-13-03','TSF-13-05','TSF-13-06')

data <- list()
for (i in 1:length(samples)) {
  print(samples[i])
  data[[i]] = read.table(paste(samples[i], ".mCount.txt", sep=""), col.names = c("chr","pos","N","X")) # see DSS documentation for the format requirement
}

BSdat <- makeBSseqData(data, samples) #  CpGs
BSdat <- sort(BSdat)

### Make a design table
design <- data.frame(morph = c(rep("White",7),
                               rep("Tan",5)),
                     age = c(rep("Chick",3), rep("Adult",4),
                             rep("Chick",2), rep("Adult",3)))
design$morph <- factor(design$morph, levels = c("Tan","White"))
design$age <- factor(design$age, levels = c("Adult","Chick"))
row.names(design) <- samples
pData(BSdat) <- design

### Filter CpGs
keepLoci <- which(DelayedMatrixStats::rowSums2(getCoverage(BSdat, type="Cov") == 0) == 0)

length(keepLoci)
BSdat.filtered <- BSdat[keepLoci,]

##### Identify DMLs #####
DMLfit = DMLfit.multiFactor(BSdat.filtered, design=design, formula=~morph+age)
coefs <- colnames(DMLfit$X)
statL <- list() # list of stats
pvalsL <- list() # list of pvalues
fdrsL <- list() # list of qvalues
bonfL <- list()
for (coef in coefs[2:length(coefs)]) {
  print(coef)
  DMLtest = DMLtest.multiFactor(DMLfit, coef=coef) # DMLfit$X can check the term to test
  if (coef == "morphWhite") {
    coef = "morph"
  } else {
    coef = "age"
  }
  statL[[paste0('stat_',coef)]] = -DMLtest$stat
  pvalsL[[paste0('pval_',coef)]] = DMLtest$pvals
  fdrsL[[paste0('fdr_',coef)]] = DMLtest$fdrs
  bonfL[[paste0('bonf_',coef)]] = p.adjust(DMLtest$pvals, method = "bonferroni")
  bonf <- bonfL[[paste0('bonf_',coef)]]
}

##### Output all CpGs #####
methy.df <- data.frame(DMLtest)
methy.df$start <- methy.df$pos
methy.df$end <- methy.df$pos
methy.df <- subset(methy.df, select=c("chr", "start", "end"))

Mean <- rownames(design)
Tan <- rownames(subset(design, morph == "Tan"))
White <- rownames(subset(design, morph == "White"))
Adult <- rownames(subset(design, age == "Adult"))
Chick <- rownames(subset(design, age == "Chick"))
gps <- c("Mean","Tan","White","Adult","Chick")


for (gp in gps) {
  print(gp)
  print(get(gp))
  methy.df[,gp] <- apply(getMeth(BSdat.filtered, region=methy.df, type = "raw", what = "perRegion")[,get(gp)], FUN=mean, MARGIN=1,na.rm=TRUE)
}

methy.df <- cbind(methy.df, 
                  diff_morph = methy.df$Tan - methy.df$White,
                  diff_age = methy.df$Adult - methy.df$Chick)
methyDat <- getMeth(BSdat.filtered, region=methy.df, type = "raw", what = "perRegion")[,samples]
methy.df <- cbind(methy.df, as.data.frame(statL), as.data.frame(pvalsL), as.data.frame(fdrsL), as.data.frame(bonfL), methyDat)

write.table(methy.df, "methyDSS.txt", row.names = F, quote = F, sep = "\t")
saveRDS(methy.df, "methyDSS.rds")
save.image("dss_morph.RData")

