## Filenames from JA's and VC's experiments

fnames.ja230 <- scan("JA230_list.txt", "a")
fnames.ja230.int <- scan("JA230_list_randomly_interleaved.txt", "a")
fnames.ja187ex <- scan("JA187_extras_list.txt", "a")
fnames.vc <- scan("VC_extended_list.txt", "a")
fnames.ja230.not.int <- fnames.ja230[is.na(match(fnames.ja230, fnames.ja230.int))]

## JA

source("ICdualsound_analysis-new.R")

#for(fname in fnames.ja230.int){
#  try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, plot = FALSE, outfile = "hmm-ja230-int-ABeqlev.txt", start = 0, end = 600, bw = 50))
#}

#for(fname in fnames.ja230.not.int){
#  try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, AB.eqlevel = TRUE, go.by.soff = FALSE, plot = FALSE, outfile = "hmm-ja230-not-int-ABeqlev.txt", start = 0, end = 600, bw = 50))
#}

#for(fname in fnames.ja187ex){
#  try(hmm.JA(fname, on.reward = TRUE, match.level = FALSE, go.by.soff = FALSE, plot = FALSE, outfile = "hmm-ja187ex.txt", start = 0, end = 600, bw = 50))
#}

#for(fname in fnames.vc){
  #try(hmm.VC(fname, 1, on.reward = TRUE, match.level = FALSE, go.by.soff = FALSE, plot = FALSE, outfile = "hmm-vc-end600.txt", start = 0, end = 600, bw = 50))
  #try(hmm.VC(fname, 2, on.reward = TRUE, match.level = FALSE, go.by.soff = FALSE, plot = FALSE, outfile = "hmm-vc-end600.txt", start = 0, end = 600, bw = 50))
  #try(hmm.VC(fname, 1, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, plot = FALSE, outfile = "hmm-vc.txt", start = 0, end = 600, bw = 50))
  #try(hmm.VC(fname, 2, on.reward = TRUE, match.level = FALSE, go.by.soff = TRUE, plot = FALSE, outfile = "hmm-vc.txt", start = 0, end = 600, bw = 50))
#}


#pdf(height = 8, width = 5, file = "poi-ja230-int-ABeqlev.pdf")
#for(fname in fnames.ja230.int){
#  try(poi.JA(fname, TRUE, FALSE, TRUE, outfile = "bfs-ja230-int-ABeqlev.txt", 0, 600))
#}
#dev.off()

#pdf(height = 8, width = 5, file = "poi-ja230-not-int-ABeqlev.pdf")
#for(fname in fnames.ja230.not.int){
#  try(poi.JA(fname, TRUE, FALSE, TRUE, outfile = "bfs-ja230-not-int-ABeqlev.txt", 0, 600))
#}
#dev.off()

#pdf(height = 8, width = 5, file = "poi-ja187ex.pdf")
#for(fname in fnames.ja187ex){
#  try(poi.JA(fname, TRUE, FALSE, outfile = "bfs-ja187ex.txt", 0, 600))
#}
#dev.off()


pdf(height = 8, width = 5, file = "poi-vc-end600.pdf")
for(fname in fnames.vc){
  try(poi.VC(fname, cell=1, on.reward=TRUE, match.level=TRUE, outfile="bfs-vc-end600.txt", start=0, end=600, go.by.soff=FALSE))
  try(poi.VC(fname, cell=2, on.reward=TRUE, match.level=TRUE, outfile="bfs-vc-end600.txt", start=0, end=600, go.by.soff=FALSE))
  #try(poi.VC(fname, 1, TRUE, TRUE, "bfs-vc.txt", 0, 600, TRUE))
  #try(poi.VC(fname, 2, TRUE, TRUE, "bfs-vc.txt", 0, 600, TRUE))
}
dev.off()


q("no")
