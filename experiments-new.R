## Filenames from JA's and VC's experiments

fnames.ja230 <- scan("JA230_list.txt", "a")
fnames.ja230.int <- scan("JA230_list_randomly_interleaved.txt", "a")
fnames.ja187ex <- scan("JA187_extras_list.txt", "a")
fnames.vc <- scan("VC_extended_list.txt", "a")
fnames.ja230.not.int <- fnames.ja230[is.na(match(fnames.ja230, fnames.ja230.int))]

## Codes needed
source("ICdualsound_analysis-new.R")
require(parallel)

## Real and synthetic control duals
dual.types <- c("realduals", "faux.allA", "faux.allB", "faux.switchingAB", "faux.averagedAB")

### VC recordings

fn.hmm.vc <- function(ty, fname, site) {
  outf <- paste("hmm-new-vc-", ty, ".txt", sep = "")
  switch(ty,
         realduals = try(hmm.VC(fname, site, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.allA = try(hmm.VC(fname, site, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, faux.dual.mix = TRUE, faux.alpha = 1, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.allB = try(hmm.VC(fname, site, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, faux.dual.mix = TRUE, faux.alpha = 0, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.switchingAB = try(hmm.VC(fname, site, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, faux.dual.swi = TRUE, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.averagedAB = try(hmm.VC(fname, site, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, faux.dual.int = TRUE, outfile = outf, start = 0, end = 600, bw = 25))
  )                                 
}

for(fname in fnames.vc) {
  for(site in c(1,2)){
    mclapply(dual.types, fn.hmm.vc, fname = fname, site = site, mc.cores = 1)
  }
}

### JA Recordings

fn.hmm.ja <- function(ty, fname, JA.set) {
  outf <- paste("hmm-new-", JA.set, "-ABeqlev-", ty, ".txt", sep = "")
  switch(ty,
         realduals = try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.allA = try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, faux.dual.mix = TRUE, faux.alpha = 1, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.allB = try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, faux.dual.mix = TRUE, faux.alpha = 0, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.switchingAB = try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, faux.dual.swi = TRUE, outfile = outf, start = 0, end = 600, bw = 25)),
         faux.averagedAB = try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, faux.dual.int = TRUE, outfile = outf, start = 0, end = 600, bw = 25))
  )                                 
}

for(JA.set in c("ja230.int", "ja230.not.int", "ja.187ex")){
  for(fname in fnames.ja230.int) {
    mclapply(dual.types, fn.hmm.ja, fname = fname, JA.set = JA.set, mc.cores = 5)
  }
}

q("no")

