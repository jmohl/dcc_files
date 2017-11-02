setwd("C:/Users/jtm47/Documents/MATLAB/RCode/reference")

## VC ##

bfs <- read.table("bfs-vc-end600.txt")
names(bfs) <- c("CellId", "cell", "frq", "pos", "sep", "mix", "ave", "out", "dom", "pval1", "pval2", "pvalD")
summary(bfs)
attach(bfs)

all.bf <- data.frame(Mixture = mix, Average = ave, Outside = out, Single = dom)
p.post <- exp(all.bf) / rowSums(exp(all.bf))
summary(p.post)

WinModel <- dimnames(p.post)[[2]][apply(p.post, 1, which.max)]
post.probs <- as.data.frame(p.post)
names(post.probs) <- c("PrMix", "PrAve", "PrOut", "PrSing")
post.probs$WinModels <- WinModel
post.probs$WinPr <- apply(post.probs[,1:4], 1, max)

dd <- bfs[,c(1:5,10:11)]
names(dd) <- c("CellId", "PairId", "AltFreq", "AltPos", "SepBF", "Pval1", "Pval2")
full.dd <- cbind(dd, post.probs)
write.csv(full.dd, file = "poi-vc.csv", row.names = FALSE)
detach(bfs)

poi.vc <- read.csv("poi-vc.csv")
hmm.vc <- read.table("hmm-vc-end600.txt")
names(hmm.vc) <- c("CellId", "PairId", "AltFreq", "AltPos", "MlpxRate", "PrState1", "Switch12", "Switch21", "ProbMlpx", "PvalSep", "PrBinSep", "MinSampSize", paste("Pswitch", seq(0,90,10), "-", seq(10,100,10), sep = ""))
results.vc <- merge(poi.vc, hmm.vc, all = TRUE)

results.vc$'Rule-5-5-5' <- results.vc$MinSampSize >= 5 & rowSums(results.vc[,grep("Pswitch", names(results.vc))]) >= 5
write.csv(results.vc, file = "combined-vc.csv", row.names = FALSE)


## JA-230-int ##

#bfs <- read.table("bfs-ja230-int.txt")
bfs <- read.table("bfs-ja230-int-ABeqlev.txt")
names(bfs) <- c("CellId", "frq", "pos", "sep", "mix", "ave", "out", "dom", "pval1", "pval2", "pvalD")
summary(bfs)
attach(bfs)

all.bf <- data.frame(Mixture = mix, Average = ave, Outside = out, Single = dom)
p.post <- exp(all.bf) / rowSums(exp(all.bf))
summary(p.post)

WinModel <- dimnames(p.post)[[2]][apply(p.post, 1, which.max)]
post.probs <- as.data.frame(p.post)
names(post.probs) <- c("PrMix", "PrAve", "PrOut", "PrSing")
post.probs$WinModels <- WinModel
post.probs$WinPr <- apply(post.probs[,1:4], 1, max)

dd <- bfs[,c(1:4,9:10)]
names(dd) <- c("CellId", "AltFreq", "AltPos", "SepBF", "Pval1", "Pval2")
full.dd <- cbind(dd, post.probs)
#write.csv(full.dd, file = "poi-ja230-int.csv", row.names = FALSE)
write.csv(full.dd, file = "poi-ja230-int-ABeqlev.csv", row.names = FALSE)
detach(bfs)

poi.ja <- read.csv("poi-ja230-int-ABeqlev.csv")
hmm.ja <- read.table("hmm-ja230-int-ABeqlev.txt")
names(hmm.ja) <- c("CellId", "AltFreq", "AltPos", "MlpxRate", "PrState1", "Switch12", "Switch21", "ProbMlpx", "PvalSep", "PrBinSep", "MinSampSize", paste("Pswitch", seq(0,90,10), "-", seq(10,100,10), sep = ""))
results.ja <- merge(poi.ja, hmm.ja, all = TRUE)

results.ja$'Rule-5-5-5' <- results.ja$MinSampSize >= 5 & rowSums(results.ja[,grep("Pswitch", names(results.ja))]) >= 5
write.csv(results.ja, file = "combined-ja230-int-ABeqlev.csv", row.names = FALSE)


## JA-230-not-int ##

bfs <- read.table("bfs-ja230.txt")
names(bfs) <- c("CellId", "frq", "pos", "sep", "mix", "ave", "out", "dom", "pval1", "pval2", "pvalD")
summary(bfs)
attach(bfs)

all.bf <- data.frame(Mixture = mix, Average = ave, Outside = out, Single = dom)
p.post <- exp(all.bf) / rowSums(exp(all.bf))
summary(p.post)

WinModel <- dimnames(p.post)[[2]][apply(p.post, 1, which.max)]
post.probs <- as.data.frame(p.post)
names(post.probs) <- c("PrMix", "PrAve", "PrOut", "PrSing")
post.probs$WinModels <- WinModel
post.probs$WinPr <- apply(post.probs[,1:4], 1, max)

dd <- bfs[,c(1:4,9:10)]
names(dd) <- c("CellId", "AltFreq", "AltPos", "SepBF", "Pval1", "Pval2")
full.dd <- cbind(dd, post.probs)
write.csv(full.dd, file = "poi-ja230.csv", row.names = FALSE)
detach(bfs)

poi.ja <- read.csv("poi-ja230.csv")
hmm.ja <- read.table("hmm-ja230.txt")
names(hmm.ja) <- c("CellId", "AltFreq", "AltPos", "MlpxRate", "PrState1", "Switch12", "Switch21", "ProbMlpx", "PvalSep", "PrBinSep", "MinSampSize", paste("Pswitch", seq(0,90,10), "-", seq(10,100,10), sep = ""))
results.ja <- merge(poi.ja, hmm.ja, all = TRUE)

results.ja$'Rule-5-5-5' <- results.ja$MinSampSize >= 5 & rowSums(results.ja[,grep("Pswitch", names(results.ja))]) >= 5
write.csv(results.ja, file = "combined-ja230-not-int-ABeqlev.csv", row.names = FALSE)

## JA-187ex ##

bfs <- read.table("bfs-ja187ex.txt")
names(bfs) <- c("CellId", "frq", "pos", "sep", "mix", "ave", "out", "dom", "pval1", "pval2", "pvalD")
summary(bfs)
attach(bfs)

all.bf <- data.frame(Mixture = mix, Average = ave, Outside = out, Single = dom)
p.post <- exp(all.bf) / rowSums(exp(all.bf))
summary(p.post)

WinModel <- dimnames(p.post)[[2]][apply(p.post, 1, which.max)]
post.probs <- as.data.frame(p.post)
names(post.probs) <- c("PrMix", "PrAve", "PrOut", "PrSing")
post.probs$WinModels <- WinModel
post.probs$WinPr <- apply(post.probs[,1:4], 1, max)

dd <- bfs[,c(1:4,9:10)]
names(dd) <- c("CellId", "AltFreq", "AltPos", "SepBF", "Pval1", "Pval2")
full.dd <- cbind(dd, post.probs)
write.csv(full.dd, file = "poi-ja187ex.csv", row.names = FALSE)
detach(bfs)

poi.ja <- read.csv("poi-ja187ex.csv")
hmm.ja <- read.table("hmm-ja187ex.txt")
names(hmm.ja) <- c("CellId", "AltFreq", "AltPos", "MlpxRate", "PrState1", "Switch12", "Switch21", "ProbMlpx", "PvalSep", "PrBinSep", "MinSampSize", paste("Pswitch", seq(0,90,10), "-", seq(10,100,10), sep = ""))
results.ja <- merge(poi.ja, hmm.ja, all = TRUE)

results.ja$'Rule-5-5-5' <- results.ja$MinSampSize >= 5 & rowSums(results.ja[,grep("Pswitch", names(results.ja))]) >= 5
write.csv(results.ja, file = "combined-ja187ex.csv", row.names = FALSE)


### ALL Combined ###

res1 <- read.csv("combined-vc.csv")
res2 <- read.csv("combined-ja230-int.csv")
res3 <- read.csv("combined-ja230-not-int.csv")
res4 <- read.csv("combined-ja187ex.csv")

res12 <- merge(res1, res3, by = names(res3), all = TRUE)
res123 <- merge(res12, res3, by = names(res3), all = TRUE)
res <- merge(res123, res4, by = names(res4), all = TRUE)
write.csv(res, file = "combined-all.csv", row.names = FALSE)


### ALL Interleaved ###

res1 <- read.csv("combined-vc.csv")
res2 <- read.csv("combined-ja230-int.csv")
res3 <- read.csv("combined-ja187ex.csv")

res12 <- merge(res1, res2, by = names(res2), all = TRUE)
res <- merge(res12, res3, by = names(res3), all = TRUE)
write.csv(res, file = "combined-all-int.csv", row.names = FALSE)