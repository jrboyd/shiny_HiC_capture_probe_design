source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)

full_gr = readRDS(bfcrpath(enz_bfc, "hg38,HindIII"))
# gr = sample(full_gr, 200)
gr = full_gr
seq = getSeq(Hsapiens, names = gr)
freq =  as.data.table(alphabetFrequency(seq, baseOnly=TRUE))
freq$id = gr$name
freq[, cg := (C + G) / (A + C + G + T)]
freq[, w := A + C + G + T]

full_gr$cg = freq$cg


plot(freq$cg[1:2000])
