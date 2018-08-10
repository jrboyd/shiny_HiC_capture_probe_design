library(BiocFileCache)
library(GenomicRanges)
library(data.table)
library(pbapply)
library(Gviz)

select_probes = function(probes_f, fragments_f, feature_gr, out_dir = getwd(), out_root = "selected_probes"){
    dir.create(out_dir, showWarnings = FALSE)
    # probes_dt = rbind(fread("data_PG_hist_capture/Ghule_Targets_Probes.txt"))
    probes_dt = rbindlist(lapply(probes_f, fread))
    probes_gr = GRanges(probes_dt$Coordinates)
    probes_gr$probe_id = probes_dt$ProbeID
    
    # frags_gr = rtracklayer::import.bed("data_PG_hist_capture/frags-2.bed")
    frags_gr = rtracklayer::import.bed(fragments_f)
    strand(frags_gr) = "*"
    
    # ref_bfc = BiocFileCache(cache = ".cache_ref")
    # ref_gr = readRDS(bfcrpath(ref_bfc, "hg38,PG_hist"))
    # pr_gr = promoters(ref_gr, 1000, 1000)
    
    olaps = findOverlaps(query = probes_gr, subject = frags_gr, type = "within")
    probes_inFrags_gr = probes_gr[queryHits(olaps)]
    
    sapply(width(frags_gr), function(w)min(w, 330))
    small_frags_gr = frags_gr[width(frags_gr) <= 330]
    wide_frags_gr = frags_gr[width(frags_gr) > 330]
    
    left_frags_gr = resize(wide_frags_gr, 330, fix = "start")
    right_frags_gr = resize(wide_frags_gr, 330, fix = "end")
    
    # combine capture regions for plotting
    cap_gr = reduce(c(small_frags_gr, left_frags_gr, right_frags_gr))
    
    # pick leftmost probe in left, rightmost in right, and left and right-most in small
    olaps = findOverlaps(query = probes_inFrags_gr, subject = left_frags_gr, type = "within")
    left_probes_inFrags_gr = probes_inFrags_gr[queryHits(olaps)]
    left_probes_inFrags_gr$subjectHits = subjectHits(olaps)
    
    olaps = findOverlaps(query = probes_inFrags_gr, subject = right_frags_gr, type = "within")
    right_probes_inFrags_gr = probes_inFrags_gr[queryHits(olaps)]
    right_probes_inFrags_gr$subjectHits = subjectHits(olaps)
    
    olaps = findOverlaps(query = probes_inFrags_gr, subject = small_frags_gr, type = "within")
    small_probes_inFrags_gr = probes_inFrags_gr[queryHits(olaps)]
    small_probes_inFrags_gr$subjectHits = subjectHits(olaps)
    
    #leftmost left
    left_dt = unique(as.data.table(left_probes_inFrags_gr))
    table(left_dt[, .N, by = subjectHits]$N)
    left_dt[, isLeft := start == min(start), by = .(subjectHits)]
    table(unique(left_dt)[isLeft == TRUE, .N, by = subjectHits]$N)
    leftmost_probes_inFrags_gr = GRanges(left_dt[isLeft == TRUE,])
    #rightmost right
    right_dt = unique(as.data.table(right_probes_inFrags_gr))
    table(right_dt[, .N, by = subjectHits]$N)
    right_dt[, isRight := end == max(end), by = .(subjectHits)]
    table(unique(right_dt)[isRight == TRUE, .N, by = subjectHits]$N)
    rightmost_probes_inFrags_gr = GRanges(right_dt[isRight == TRUE,])
    #leftmost small
    left_dt = unique(as.data.table(small_probes_inFrags_gr))
    table(left_dt[, .N, by = subjectHits]$N)
    left_dt[, isLeft := start == min(start), by = .(subjectHits)]
    table(unique(left_dt)[isLeft == TRUE, .N, by = subjectHits]$N)
    sm_leftmost_probes_inFrags_gr = GRanges(left_dt[isLeft == TRUE,])
    #rightmost small
    right_dt = unique(as.data.table(small_probes_inFrags_gr))
    table(right_dt[, .N, by = subjectHits]$N)
    right_dt[, isRight := end == max(end), by = .(subjectHits)]
    table(unique(right_dt)[isRight == TRUE, .N, by = subjectHits]$N)
    sm_rightmost_probes_inFrags_gr = GRanges(right_dt[isRight == TRUE,])
    
    mcols(leftmost_probes_inFrags_gr) = NULL
    mcols(rightmost_probes_inFrags_gr) = NULL
    mcols(sm_leftmost_probes_inFrags_gr) = NULL
    mcols(sm_rightmost_probes_inFrags_gr) = NULL
    
    nearestends_probes_inFrags_gr = unique(c(leftmost_probes_inFrags_gr, 
                                             rightmost_probes_inFrags_gr, 
                                             sm_leftmost_probes_inFrags_gr, 
                                             sm_rightmost_probes_inFrags_gr))
    
    # probes outside capture regions should be removed
    uncap_frags_gr = setdiff(wide_frags_gr, cap_gr)
    
    olaps = findOverlaps(query = probes_inFrags_gr, subject = cap_gr, type = "within")
    cap_probes_inFrags_gr = probes_inFrags_gr[queryHits(olaps)]
    
    cap_gr$color = "blue"
    uncap_frags_gr$color = "gray"
    
    annot_frags_gr = c(cap_gr, uncap_frags_gr)
    
    myProbeTracks_i = function(i, hgr){
        qgr = resize(frags_gr[i], width = 10000, fix = "center")
        ### setup tracks
        myProbeTracks(qgr, hgr)
    }
    
    myProbeTracks = function(qgr, hgr){
        probes_all_Track <- AnnotationTrack(range = subsetByOverlaps(probes_gr, qgr),
                                            genome = "hg38", name = "All Probes")
        probes_inFrags_Track <- AnnotationTrack(range = subsetByOverlaps(probes_inFrags_gr, qgr),
                                                genome = "hg38", name = "Probes in Fragments")
        
        probes_cap_Track = AnnotationTrack(range = subsetByOverlaps(cap_probes_inFrags_gr, qgr),
                                           genome = "hg38", name = "Probes within 330bp of Ends of Fragments")
        probes_nearest_Track = AnnotationTrack(range = subsetByOverlaps(nearestends_probes_inFrags_gr, qgr),
                                               genome = "hg38", name = "Probes Nearest Ends of Fragments")
        
        plots = list(probes_all_Track, probes_inFrags_Track, probes_cap_Track, probes_nearest_Track)
        
        if(is.null(hgr$color)) hgr$color = "black"
        ht = HighlightTrack(plots, range = subsetByOverlaps(hgr, qgr), col = "black", fill = subsetByOverlaps(hgr, qgr)$color, alpha = .3)
        plotTracks(ht, from = start(qgr), to = end(qgr), main = as.character(qgr))    
    }
    
    ### plot tracks
    todo = sample(seq_along(frags_gr), 10)
    tracks_pdf = file.path(out_dir, paste0(out_root, "_Tracks.pdf"))
    pdf(tracks_pdf)
    sapply(todo, function(i){
        print(i)
        myProbeTracks_i(i, annot_frags_gr)
    })
    dev.off()
    
    final_probes_gr = nearestends_probes_inFrags_gr
    final_sizeMb = sum(width(final_probes_gr))/10^6
    olaps = findOverlaps(final_probes_gr, frags_gr, type = "within")
    odt = as.data.table(olaps)
    nolaps = odt[, .N, by = .(subjectHits)]
    frags_gr$n_probes = 0
    frags_gr[nolaps$subjectHits]$n_probes = nolaps$N
    table(frags_gr$n_probes)
    frags_gr$score = NULL
    ref_gr
    distToFrag = as.data.table(distanceToNearest(x = feature_gr, subject = frags_gr))
    distToFrag$distance
    
    ref_dt = as.data.table(ref_gr)
    
    ref_gr$distance_to_fragment = Inf
    ref_gr[distToFrag$queryHits]$distance_to_fragment = distToFrag$distance
    
    
    
    ref_gr$frag_id = "no_fragment"
    ref_gr$n_probes = 0L
    ref_gr$frag_id = frags_gr[distToFrag$subjectHits]$name
    ref_gr$n_probes = frags_gr[distToFrag$subjectHits]$n_probes
    
    subset(ref_gr, distance_to_fragment > 10000)
    ref_gr[ref_gr$distance_to_fragment > 10000]$frag_id = "over_10kb"
    ref_gr[ref_gr$distance_to_fragment > 10000]$n_probes = 0
    
    ref_df = as.data.frame(ref_gr)
    
    bedFront = c(1:3, 13, 8, 5)
    ref_df = ref_df[, c(bedFront, setdiff(seq_len(ncol(ref_df)), bedFront))]
    
    ref_txt = file.path(out_dir, paste0(out_root, "_Reference.txt"))
    write.table(ref_df, ref_txt, sep = "\t", row.names = FALSE, quote = FALSE)
    myProbeTracks(resize(feature_gr[3], width = 10000, fix = "center"), annot_frags_gr)
    
    olaps = findOverlaps(final_probes_gr, probes_gr)
    final_probes_gr$probe_id = ""
    final_probes_gr[queryHits(olaps)]$probe_id = probes_gr[subjectHits(olaps)]$probe_id
    # subset(final_probes_gr, probe_id == "BA_j72483_000054")
    # subset(probes_gr, probe_id == "BA_j72483_000054")
    names(final_probes_gr) = final_probes_gr$probe_id
    
    probes_bed = file.path(out_dir, paste0(out_root, "_SelectedProbes.bed"))
    rtracklayer::export.bed(final_probes_gr, probes_bed)        
    
    {
        message("Final size is ", final_sizeMb, " Mb")
        message("Final probe count ", length(final_probes_gr))
        message("Probes hit ", length(subsetByOverlaps(frags_gr, final_probes_gr)), " fragments of input ", length(frags_gr))
        message("n_probes by transcript")
        print(table(ref_gr$n_probes))
        
        ref_dt = as.data.table(ref_gr)
        message("max n_probes by gene")
        print(table(ref_dt[, max(n_probes), .(gene_name)]$V1))
    }
}

bfc_ref = BiocFileCache(".cache_ref/")
ref_gr = readRDS(bfcrpath(bfc_ref, "hg38,PG_hist_round2"))
pr_gr = promoters(ref_gr, 5000, 5000)
pr_gr = resize(pr_gr, width(pr_gr) + 20000, fix = "center")

select_probes("data_PG_hist_capture_round2/Capture_Frags_rd2_Probes.txt", 
              "data_PG_hist_capture_round2/capture frags-Round2.bed", 
              pr_gr, 
              out_dir = "data_PG_hist_capture_round2/", 
              out_root = "PG_round2_")
