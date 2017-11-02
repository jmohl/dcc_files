outf <- list(VC = c("hmm-new-vc-faux.allA", "hmm-new-vc-faux.allB", 
                    "hmm-new-vc-faux.switchingAB", "hmm-new-vc-faux.averagedAB", 
                    "hmm-new-vc-realduals"),
             JA = c("hmm-new-ja230.int-ABeqlev-faux.allA", "hmm-new-ja230.int-ABeqlev-faux.allB", 
                    "hmm-new-ja230.int-ABeqlev-faux.switchingAB", "hmm-new-ja230.int-ABeqlev-faux.averagedAB", 
                    "hmm-new-ja230.int-ABeqlev-realduals",
                    "hmm-new-ja230.not.int-ABeqlev-faux.allA", "hmm-new-ja230.not.int-ABeqlev-faux.allB", 
                    "hmm-new-ja230.not.int-ABeqlev-faux.switchingAB", "hmm-new-ja230.not.int-ABeqlev-faux.averagedAB", 
                    "hmm-new-ja230.not.int-ABeqlev-realduals"),
             JM = c("hmm-new-synth_data_test-JM-faux.allA","hmm-new-synth_data_test-JM-faux.allB","hmm-new-synth_data_test-JM-faux.switchingAB",
                    "hmm-new-synth_data_test-JM-faux.averagedAB","hmm-new-synth_data_test-JM-realduals")
)

for(fname in outf[["VC"]]){
   ff <- read.table(paste(fname, ".txt", sep = ""))
   ffname <- c("Cell", "Site", "AltFrq", "AltPos",  
               "RateSwi", "RateDiv", "RateStartA", "RateA2B", "RateB2A",  
               "SwiMajor", "DivMajor",  "dupli1", "dupli2", "redundant1",
               "PvalSepRankTest",  "PvalSepPoiTest",  "MinSampleSingle",  
               "redundant2", "redundant3",  "Swi40.60",  "Swi60.80",  "Swi80.100",  
               "redundant4", "redundant5",  "Div40.60",  "Div60.80",  "Div80.100",
               "redundant6", "redundant7", "redundant8", "redundant9", "redundant10")         
   
   names(ff) <- ffname
   ff <- ff[, -unique(c(grep("redundant", ffname), grep("dupli", ffname)))]
   write.csv(ff, file = paste(fname, ".csv", sep = ""), row.names = FALSE)
}

for(fname in outf[["JA"]]){
   ff <- read.table(paste(fname, ".txt", sep = ""))
   ffname <- c("Cell", "AltFrq", "AltPos",  
               "RateSwi", "RateDiv", "RateStartA", "RateA2B", "RateB2A",  
               "SwiMajor", "DivMajor",  "dupli1", "dupli2", "redundant1",
               "PvalSepRankTest",  "PvalSepPoiTest",  "MinSampleSingle",  
               "redundant2", "redundant3",  "Swi40.60",  "Swi60.80",  "Swi80.100",  
               "redundant4", "redundant5",  "Div40.60",  "Div60.80",  "Div80.100",
               "redundant6", "redundant7", "redundant8", "redundant9", "redundant10")         
   
   names(ff) <- ffname
   ff <- ff[, -unique(c(grep("redundant", ffname), grep("dupli", ffname)))]
   write.csv(ff, file = paste(fname, ".csv", sep = ""), row.names = FALSE)
}

for(fname in outf[["JM"]]){
   ff <- read.table(paste(fname, ".txt", sep = ""))
   ffname <- c("Cell", "AltFrq", "AltPos",  
               "RateSwi", "RateDiv", "RateStartA", "RateA2B", "RateB2A",  
               "SwiMajor", "DivMajor",  "dupli1", "dupli2", "redundant1",
               "PvalSepRankTest",  "PvalSepPoiTest",  "MinSampleSingle",  
               "redundant2", "redundant3",  "Swi40.60",  "Swi60.80",  "Swi80.100",  
               "redundant4", "redundant5",  "Div40.60",  "Div60.80",  "Div80.100",
               "redundant6", "redundant7", "redundant8", "redundant9", "redundant10")         
   
   names(ff) <- ffname
   ff <- ff[, -unique(c(grep("redundant", ffname), grep("dupli", ffname)))]
   write.csv(ff, file = paste(fname, ".csv", sep = ""), row.names = FALSE)
}

##BFS processor

bfs <- read.table("bfs-synth_data50v20_weak_mplx-JM.txt")
names(bfs) <- c("CellId","frq", "pos", "sep", "mix", "ave", "out", "dom", "pval1", "pval2", "pvalD")

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
names(dd) <- c("CellId", "AltFreq","AltPos", "SepBF", "Pval1", "Pval2")
full.dd <- cbind(dd, post.probs)
write.csv(full.dd, file = "poi-synth_data50v20_mplx_range-jm.csv", row.names = FALSE)
detach(bfs)