check_ref_bfc = function(rnam){
    nrow(bfcquery(ref_bfc, rnam, field="rname")) > 0
}

check_enz_bfc = function(rnam){
    nrow(bfcquery(enz_bfc, rnam, field="rname")) > 0
}

bfcrcheck = function(x, rnam){
    nrow(bfcquery(x, rnam, field="rname")) > 0
}

bfcrget = function(x, rnam){
    q = bfcquery(x, paste0("^", rnam, "$"), field="rname")
    if(nrow(q) < 1){
        q = bfcquery(x, rnam, field="rname")
    }
    if(nrow(q) > 1){
        stop("too many results")
    }
    if(nrow(q) < 1){
        stop("no matches")
    }
    paste0(bfccache(x), "/", q$rpath)
}