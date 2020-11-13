# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg38")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(data.table)

full_gr = readRDS(bfcrpath(enz_bfc, "hg38,HindIII"))


enz_all = bfcinfo(enz_bfc)
for(rname in enz_all$rname){
    f = bfcrpath(enz_bfc, rname)
    message(rname, " ", f)
    full_gr = readRDS(f)
    if(!"GRanges" %in% class(full_gr)){
        message("skip, not GRanges")
        next
    }
    if(is.null(full_gr$CG)){
        message("need CG")
        gr = full_gr
        
        gen = strsplit(rname, ",")[[1]][1]
        if(gen == "hg38"){
            seq = getSeq(Hsapiens, names = gr)    
        }else if(gen == "mm10"){
            seq = getSeq(Mmusculus, names = gr)
        }else{
            stop("bad gen")
        }
        
        freq =  as.data.table(alphabetFrequency(seq, baseOnly=TRUE))
        freq$id = gr$name
        freq[, cg := (C + G) / (A + C + G + T)]
        freq[, w := A + C + G + T]
        
        full_gr$CG = freq$cg
        saveRDS(full_gr, f)
    }else{
        message("has CG")
    }
}

# gr = sample(full_gr, 200)



plot(freq$cg[1:2000])
