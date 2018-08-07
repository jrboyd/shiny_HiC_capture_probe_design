ucsc_probe_tracks = function(enzyme, 
                             enz_full, 
                             enz_filtered, 
                             enz_selected, 
                             annotation_full, 
                             annotation_filtered, 
                             promoters_filtered,
                             promoters_missed,
                             out_dir = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/joeboyd_files/probe_app/",
                             cfg_f = paste0(out_dir, "track_config_probe_design.txt"),
                             uid = round(runif(1)*10^6)
){
    file.remove(dir(out_dir, full.names = TRUE))
    writeTrack = function(gr, desc, bedf, bbf = sub("\\.bed$", "\\.bb", bedf), color = "black"){
        colrgb = paste(col2rgb(color)[,1], collapse = ",")
        rtracklayer::export.bed(gr, bedf)
        system(paste("bedSort", bedf, paste0(bedf, ".sorted")), intern = TRUE)
        system(paste("bedToBigBed", paste0(bedf, ".sorted"), "~/hg38_chrsizes.txt", bbf), intern = TRUE)
        makeTrack(file = bbf, type = "bigBed", name = desc, description = desc, color = colrgb)
    }
    
    annotation_full$score = 0
    annotation_full$name = annotation_full$gene_name
    annotation_full$strand = "*"
    annotation_full = sortSeqlevels(GRanges(annotation_full))
    annotation_full = sort(annotation_full)
    
    annotation_filtered$score = 0
    annotation_filtered$name = annotation_filtered$gene_name
    annotation_filtered$strand = "*"
    annotation_filtered = sortSeqlevels(GRanges(annotation_filtered))
    annotation_filtered = sort(annotation_filtered)
    
    promoters_filtered = as.data.frame(promoters_filtered)
    promoters_filtered$score = 0
    promoters_filtered$name = promoters_filtered$gene_name
    promoters_filtered$strand = "*"
    promoters_filtered = sortSeqlevels(GRanges(promoters_filtered))
    promoters_filtered = sort(promoters_filtered)
    
    promoters_missed = as.data.frame(promoters_missed)
    if(nrow(promoters_missed) > 0){
        promoters_missed$score = 0    
        promoters_missed$strand = "*"
        promoters_missed$name = promoters_missed$gene_name
        promoters_missed = sortSeqlevels(GRanges(promoters_missed))
        promoters_missed = sort(promoters_missed)
    }else{
        promoters_missed = GRanges("chr21", IRanges(1, 2))
    }
    
    
    
    track_lines = c(
        writeTrack(enz_full, paste(enzyme, "Full"), paste0(out_dir, uid, "_", enzyme, "_enz_full.bed"), color = "gray"),
        writeTrack(enz_filtered, paste(enzyme, "Filtered"), paste0(out_dir, uid, "_", enzyme, "_enz_filtered.bed"), color = "black"),
        writeTrack(enz_selected, paste(enzyme, "Selected"), paste0(out_dir, uid, "_", enzyme, "_enz_selected.bed"), color = "darkgreen"),
        writeTrack(annotation_full, paste("Transcripts Full"), paste0(out_dir, uid, "_annot_full.bed"), color = "gray"),
        writeTrack(annotation_filtered, paste("Transcripts Filtered"), paste0(out_dir, uid, "_annot_filtered.bed"), color = "black"),
        writeTrack(promoters_filtered, paste("Promoters Filtered"), paste0(out_dir, uid, "_prom_filtered.bed"), color = "black"),
        writeTrack(promoters_missed, paste("Promoters Missed"), paste0(out_dir, uid, "_", enzyme, "_prom_missed.bed"), color = "red")
    )
    
    writeLines(track_lines, cfg_f)
    return(cfg_f)
}

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
    # browser()
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
    # paste0(bfccache(x), "/", q$rpath)
    q$rpath
}

TEMPLATE = paste("track name=\"NAME\" description=\"DESCRIPTION\"", 
                 "visibility=VISIBILITY",
                 "color=COLOR bigDataUrl=URL type=TYPE")
FILE_ROOT = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks"
URL_ROOT = "https://galaxy.med.uvm.edu/static/UCSCtracks"
makeTrack = function(file, type, name, description,  
                     color, visibility = "pack"){
    url = sub(FILE_ROOT, URL_ROOT, file)
    new_track = TEMPLATE
    new_track = sub("TYPE", type, new_track)
    new_track = sub("URL", url, new_track)
    new_track = sub("NAME", name, new_track)
    new_track = sub("DESCRIPTION", description, new_track)
    new_track = sub("VISIBILITY", visibility, new_track)
    new_track = sub("COLOR", color, new_track)
    return(new_track)
}
