#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
DEV=!basename(dirname(getwd())) == "ShinyApps"
DEV = FALSE
library(shiny)
library(shinycssloaders)
library(BiocFileCache)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(cowplot)
source("functions_app.R")
ref_urls = list("hg38" = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz")
ref_files = list()
raw_bfc = BiocFileCache(cache = ".cache")
ref_bfc = BiocFileCache(cache = ".cache_ref")
enz_bfc = BiocFileCache(cache = ".cache_enzyme")

#use biocfile cache to store reference as text and pre-loaded GRanges if needed
for(i in seq_along(ref_urls)){
    fpath = bfcrpath(x = raw_bfc, ref_urls[[i]], rnames = names(ref_urls)[i]) 
    ref_files[[names(ref_urls)[i]]] = fpath
    
    rnam_full = paste(names(ref_urls)[i], "FULL", sep = ",")
    if(!bfcrcheck(ref_bfc, rnam_full)){
        ref = rtracklayer::import.gff(fpath, format = "gtf", feature.type = "transcript")
        saveRDS(ref, bfcnew(ref_bfc, rnam_full))
    }
}

if(DEV & !exists("tmp_ref")){
    tmp_ref = readRDS(bfcrget(ref_bfc, "hg38,FULL"))
    tmp_ref = subset(tmp_ref, seqnames == "chr21")
} 
enz_files = c("hg38,HindIII" = "~/HiC-Pro/data/frags_hg38canon_HindIII.bed",
              "hg38,DpnII" = "~/HiC-Pro/data/frags_hg38canon_MboI.bed")

for(i in seq_along(enz_files)){
    rname = names(enz_files)[i]
    if(!bfcrcheck(enz_bfc, rname)){
        enz = rtracklayer::import.bed(enz_files[i])
        saveRDS(enz, bfcnew(enz_bfc, rname))
    }
}

# if(DEV & !exists("tmp_enz")){
#     tmp_enz = readRDS(bfcrpath(enz_bfc, "hg38,HindIII"))
#     tmp_enz = subset(tmp_enz, seqnames == "chr21")
# } 
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("HiC capture probe design"),
    h5("Move down the page, working through each step."),
    # Sidebar with a slider input for number of bins 
    h2("Step1 - Select Reference"),
    hr(),
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(selectInput("selectGenome", "Genome", choices = rev(bfcinfo(ref_bfc)$rname)), width = 6),
                column(selectInput("selectCutter", "Enzyme", choices = c("DpnII", "HindIII", "HindIII_PG_hist_round2"), selected = "HindIII_PG_hist_round2"), width = 6)
                
            ),
            fluidRow(
                column(withSpinner(htmlOutput("htmlTxCountRaw")),
                       withSpinner(htmlOutput("htmlFragCountRaw")), width = 6)
            )
        ), mainPanel(
            
        )
    ),
    h2("Step2 - Refine Annotation"),
    hr(),
    sidebarLayout(
        sidebarPanel(
            sliderInput("sliderMinFragSize", 
                        "Minimum Fragment Sizes", 
                        min = 0, max = 1000, 
                        step = 10, 
                        value = 200),
            sliderInput("sliderMaxFragSize", 
                        "Maximum Fragment Sizes", 
                        min = 5000, max = 200000, 
                        step = 5000, 
                        value = 10000),
            sliderInput("sliderCGrange",
                        "CG content",
                        min = 0, max = 1,
                        step = 0.01,
                        value = c(0, 1)),
            fluidRow(
                column(withSpinner(htmlOutput("htmlTxCountFilter")),
                       withSpinner(htmlOutput("htmlFragCountFilter")), width = 4),
                column(withSpinner(plotOutput("plotFragSize", width = "320px", height = "200px")), width = 8)
            )
            
        ), mainPanel(
            fluidRow(
                column(width = 3,
                       actionButton("btnRestart", "Start Over"),
                       hr(),
                       actionButton("btnRemoveSelected", "Remove Selected"),
                       actionButton("btnLimitToSelected", "Limit To Selected"),
                       hr(),
                       actionButton("btnRemoveVisible", "Remove Filtered"),
                       actionButton("btnLimitToVisible", "Limit To Filtered")
                ),
                column(width = 9,
                       tags$h4("Use table to filter target gene set"),
                       tags$h6("Recommend gene_type and transcript_support_level"),
                       withSpinner(DT::dataTableOutput(outputId = "filterTable"))
                )
            )
        )
    ),
    h2("Step3 - Filter Fragments"),
    hr(),
    htmlOutput("lockStatus"),
    fluidRow(
        actionButton("btnLockRef", "Lock Reference", icon = icon("unlock")),
        actionButton("btnLockEnz", "Lock Fragments", icon = icon("unlock"))
    ),
    htmlOutput("unlockStep4"),
    br(),
    br(),
    br()
)

server <- function(input, output, session) {
    
    rvTSScount = reactiveVal()
    rvAnnotRaw = reactiveVal()
    rvAnnotFiltering = reactiveVal()
    rvAnnotLock = reactiveVal()
    rvAnnotExt = reactiveVal()
    rvAnnotMissed = reactiveVal()
    # rvAnnotFinal = reactiveVal()
    
    rvEnzRaw = reactiveVal()
    rvEnzFiltering = reactiveVal()
    rvEnzLock = reactiveVal()
    rvEnzSelected = reactiveVal()
    
    rvOlapDt = reactiveVal()
    
    observeEvent(
        eventExpr = {
            input$selectCutter
            input$selectGenome
        }, 
        handlerExpr = {
            req(input$selectCutter)
            input$selectCutter
            fpath = enz_files[input$selectCutter]
            gen = sub(",.+", "", input$selectGenome)
            rname = paste(gen, input$selectCutter, sep = ",")
            # browser()
            gr = readRDS(bfcrget(enz_bfc, rname))
            if(DEV){
                gr = subset(gr, seqnames == "chr21")    
            }
            rvEnzRaw(gr)
        }
        
    )
    
    #any change to filtering invalidates locked
    observeEvent(
        eventExpr = {
            rvAnnotFiltering()
        },
        handlerExpr = {
            rvAnnotLock(NULL)
            updateActionButton(session, "btnLockRef", icon = icon("unlock"))
        }
    )
    
    #handle lock ref press
    observeEvent(
        eventExpr = {
            input$btnLockRef
        },
        handlerExpr = {
            rvAnnotLock(rvAnnotFiltering())
            updateActionButton(session, "btnLockRef", icon = icon("lock"))
        }
    )
    
    #any change to filtering invalidates locked
    observeEvent(
        eventExpr = {
            rvEnzFiltering()
        },
        handlerExpr = {
            rvEnzLock(NULL)
            updateActionButton(session, "btnLockEnz", icon = icon("unlock"))
        }
    )
    
    #handle lock enz press
    observeEvent(
        eventExpr = {
            input$btnLockEnz
        },
        handlerExpr = {
            rvEnzLock(rvEnzFiltering())
            updateActionButton(session, "btnLockEnz", icon = icon("lock"))
        }
    )
    
    
    #if selected reference changes, update rvAnnotRaw
    observeEvent(
        eventExpr = {
            input$selectGenome
        }, 
        handlerExpr = {
            # rname = paste(input$selectGenome, "FULL", sep = ",")
            rname = input$selectGenome
            # browser()
            gr = readRDS(bfcrget(ref_bfc, rname))
            if(DEV){
                gr = subset(gr, seqnames == "chr21")
            }
            rvAnnotRaw(as.data.frame(gr))
        }
    )
    
    #if rvAnnotRaw changes or start over, init rvAnnotFiltering
    observeEvent(
        eventExpr = {
            input$btnRestart
            rvAnnotRaw()
        }, 
        handlerExpr = {
            rvAnnotFiltering(rvAnnotRaw())
        }
    )
    
    observeEvent(input$btnRemoveSelected, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        
        rsel = sort(as.integer(input$filterTable_rows_selected))
        rvis = setdiff(rvis, rsel)
        if(length(rsel) < 1){
            showNotification("Empty selection.", type = "error")
        }else{
            rvAnnotFiltering(rvAnnotFiltering()[rvis,])
            showNotification(paste("Removed", length(rsel), "entries."))
        }
    })
    
    observeEvent(input$btnLimitToSelected, {
        rsel = sort(as.integer(input$filterTable_rows_selected))
        nr = nrow(rvAnnotFiltering())
        if(length(rsel) < 1){
            showNotification("Empty selection.", type = "error")
        }else{
            rvAnnotFiltering(rvAnnotFiltering()[rsel,])
            showNotification(paste("Removed", nr - length(rsel), "entries."))
        }
    })
    
    observeEvent(input$btnRemoveVisible, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        nr = nrow(rvAnnotFiltering())
        if(length(rvis) < 1){
            showNotification("Would be empty.", type = "error")
        }else if(length(rvis) == nr){
            showNotification("No filters.", type = "error")
        }else{
            rvAnnotFiltering(rvAnnotFiltering()[-rvis,])
            showNotification(paste("Removed", length(rvis), "entries."))
        }
        
    })
    
    observeEvent(input$btnLimitToVisible, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        nr = nrow(rvAnnotFiltering())
        if(length(rvis) < 1){
            showNotification("Would be empty.", type = "error")
        }else if(length(rvis) == nr){
            showNotification("No filters.", type = "error")
        }else{
            rvAnnotFiltering(rvAnnotFiltering()[rvis,])
            showNotification(paste("Removed", nr - length(rvis), "entries."))
        }
        
    })
    
    #if TSS size parameter or rvAnnotFiltering change, update rvAnnotExt
    observeEvent(
        eventExpr = 
        {
            rvAnnotLock()
            input$numericTSSdownstream
            input$numericTSSupstream
            input$numericMaxSearch
        }, 
        handlerExpr = {
            if(is.null(rvAnnotLock())){
                rvAnnotExt(rvAnnotLock())
            }else{
                gr = GRanges(rvAnnotLock())
                gr = promoters(gr, 
                               upstream = input$numericTSSupstream, 
                               downstream = input$numericTSSdownstream)
                start(gr) = start(gr) - input$numericMaxSearch
                end(gr) = end(gr) + input$numericMaxSearch
                rvAnnotExt(gr)
            }
        }
    )
    
    rvConfigFile = reactiveVal()
    
    observe({
        req(input$selectCutter)
        req(rvEnzSelected())
        req(rvAnnotMissed())
        cfg_f = ucsc_probe_tracks(
            enzyme = input$selectCutter, 
            enz_full = rvEnzRaw(), 
            enz_filtered = rvEnzFiltering(), 
            enz_selected = rvEnzSelected(), 
            annotation_full = rvAnnotRaw(), 
            annotation_filtered = rvAnnotFiltering(), 
            promoters_filtered = rvAnnotExt(),
            promoters_missed = rvAnnotMissed()
        )
        rvConfigFile(cfg_f)
    })
    
    output$lockStatus = renderUI(
        if(is.null(rvEnzLock()) || is.null(rvAnnotLock())){
            h5("You must lock in gene and fragment references to proceed!", style="color:red;")    
        }
    )
    
    output$unlockStep4 = renderUI(
        if(!(is.null(rvEnzLock()) || is.null(rvAnnotLock()))){
            div(
                h2("Step4 - Promoter and Search Parameters"),
                hr(),
                sidebarLayout(
                    sidebarPanel(
                        h5("TSS Extension"),
                        fluidRow(
                            column(numericInput("numericTSSdownstream", "Downstream", value = 5000, min = 0, max = 100000, step = 500), width = 6),
                            column(numericInput("numericTSSupstream", "Upstream", value = 5000, min = 0, max = 100000, step = 500), width = 6)
                        ),
                        numericInput("numericMaxSearch", "Max Search Extension", value = 0, min = 0, max = 100000, step = 500),
                        hr(),
                        h5("Design Parameter"),
                        sliderInput("sliderPerTSS", "Fragments per TSS", value = 1, min = 1, max = 20, step = 1),
                        sliderInput("sliderPerFragment", "Probes per fragment", value = 2, min = 1, max = 10, step = 1)
                    ), mainPanel(
                        fluidRow(
                            column(width = 3,
                                   withSpinner(plotOutput("plotTSShits", width = "240px", height = "200px")),
                                   withSpinner(plotOutput("plotOlapDist", width = "240px", height = "200px"))
                            ), 
                            column(width = 9,
                                   tags$h4("Use table to assess genes that have been missed."),
                                   tags$h6("Also, use the 'To UCSC' link to visualize"),
                                   uiOutput("ucscLink"),
                                   withSpinner(DT::dataTableOutput("missedTssTable"))
                            ))
                    )
                ),
                withSpinner(htmlOutput("calcMb")),
                downloadButton("dlFrags", label = "Download Fragments"),
                downloadButton("dlFragsPlusSeq", label = "Download Fragments With Sequence (Slow)")
                
            )
        }
    )
    
    output$htmlTxCountRaw = renderUI(
        p(paste("Transcript Count:", nrow(rvAnnotRaw())))
    )
    output$htmlFragCountRaw = renderUI(
        p(paste("Fragments Count:", length(rvEnzRaw())))
    )
    
    output$htmlTxCountFilter = renderUI(
        p(paste("Transcript Count:", nrow(rvAnnotFiltering())))
    )
    output$htmlFragCountFilter = renderUI(
        p(paste("Fragments Count:", length(rvEnzFiltering())))
    )
    
    output$filterTable = DT::renderDataTable({
        DT::datatable(rvAnnotFiltering(), filter = "top", options = list(pageLength = 25, scrollX = T)) 
    }) 
    
    output$missedTssTable = DT::renderDataTable({
        req(rvAnnotExt())
        req(rvEnzSelected())
        annot_gr = rvAnnotExt()
        enz_gr = rvEnzSelected()
        missed_gr = subsetByOverlaps(annot_gr, enz_gr, invert = TRUE, ignore.strand=TRUE)
        rvAnnotMissed(missed_gr)
        DT::datatable(as.data.frame(missed_gr), filter = "top", options = list(stateSave = TRUE, pageLength = 25, scrollX = T))
    })
    
    
    output$calcMb = renderUI({
        req(rvAnnotLock())
        req(rvEnzSelected())
        rvis = sort(as.integer(input$filterTable_rows_all))
        ntss = nrow(rvAnnotLock()[rvis,])
        total = 120 * input$sliderPerFragment * input$sliderPerTSS * ntss
        total_str = round(total / 10^6, digits = 2)
        
        nr_enz = length(rvEnzSelected())
        nr_total = nr_enz*120 * input$sliderPerFragment
        nr_total_str = round(nr_total / 10^6, digits = 2)
        tags$p(
            fluidRow(
                column(width = 3,
                       tags$h3("Probes: Maximum Size"),
                       tags$hr(),
                       tags$span(paste("Possible Total:", total_str, "Mb"), style="font-weight:bold"),
                       tags$br(),
                       tags$span(paste("TSS Count:", ntss)),
                       tags$br(),
                       tags$span(paste("Fragments Per TSS:", input$sliderPerTSS)),
                       tags$br(),
                       tags$span(paste("Total Fragments:", input$sliderPerTSS * ntss)),
                       tags$br(),
                       tags$span(paste("Probes Per Fragment:", input$sliderPerFragment)),
                       tags$br(),
                       tags$span(paste("Total Probes:", input$sliderPerFragment * input$sliderPerTSS * ntss)),
                       tags$br(),
                       tags$span("Probe Size: 120 bp"),
                       tags$br()
                       
                ),
                column(width = 3,
                       tags$h3("Probes: Estimated Size"),
                       tags$hr(),
                       tags$span(paste("Non-Redundant Total:", nr_total_str, "Mb"), style="font-weight:bold"),
                       tags$br(),
                       tags$span(paste("non-redundant Fragments:", nr_enz)),
                       tags$br(),
                       tags$span(paste("Probes Per Fragment:", input$sliderPerFragment)),
                       tags$br(),
                       tags$span("Probe Size: 120 bp"),
                       tags$br()
                       
                ),
                column(width = 3,
                       tags$h3("Regions: Actual Size"),
                       tags$hr(),
                       tags$span(paste("Fragment Size Total:", round(sum(width(rvEnzSelected()))/10^6, 1), "Mb"), style="font-weight:bold"),
                       tags$br(),
                       tags$span(paste("non-redundant Fragments:", nr_enz)),
                       tags$br(),#need average fragment size
                       tags$span(paste("average fragment size:", paste(round(mean(width(rvEnzSelected()))), "bp"))),
                       tags$br()
                       
                )
                
            )
            
        )
    })
    
    output$plotFragSize = renderPlot({
        req(rvEnzRaw())
        req(input$sliderMinFragSize)
        req(input$sliderMaxFragSize)
        req(input$sliderCGrange)
        gr = rvEnzRaw()
        
        k = width(gr) >= input$sliderMinFragSize & width(gr) <= input$sliderMaxFragSize
        gr = gr[k]
        
        k = gr$CG >= min(input$sliderCGrange) & gr$CG <= max(input$sliderCGrange)
        gr = gr[k]
        
        rvEnzFiltering(gr)
        
        gr = sample(gr, min(length(gr), 5000))
        df = data.frame(width = width(gr))
        qmax = quantile(df$width, .999)
        p = ggplot(df[df$width < qmax,, drop = FALSE], aes(x = width)) + 
            geom_histogram(bins = 100) +
            coord_cartesian(xlim = c(0, max(df$width)))
        bp = ggplot2::ggplot_build(p)
        xrng = bp$layout$panel_ranges[[1]]$x.range
        yrng = bp$layout$panel_ranges[[1]]$y.range
        p + annotate("text", 
                     x = mean(xrng), 
                     y = mean(yrng), 
                     hjust = 0, 
                     label = paste0("mean ", round(mean(df$width)), 
                                    "\nmedian ", median(df$width))) +
            labs(title = "Fragment Size Distribution") + 
            cowplot::theme_cowplot()
    })
    
    output$plotTSShits = renderPlot({
        req(rvAnnotExt())
        req(rvEnzLock())
        ref_gr = rvAnnotExt()
        mcols(ref_gr) = NULL
        ref_gr = unique(ref_gr)
        ref_gr$id = paste0("tss_", seq_along(ref_gr))
        # end(ref_gr) = end(ref_gr) + input$numericMaxSearch
        # start(ref_gr) = start(ref_gr) - input$numericMaxSearch
        enz_gr = rvEnzLock()
        
        olaps = as.data.table(findOverlaps(ref_gr, enz_gr, ignore.strand = TRUE))
        
        olap_dt = cbind(as.data.table(ref_gr[olaps$queryHits]), as.data.table(enz_gr[olaps$subjectHits]))
        nc = ncol(as.data.frame(head(ref_gr)))
        colnames(olap_dt)[1:nc] = paste0("ref_", colnames(olap_dt)[1:nc])
        
        # olap_dt[, ref_start := (ref_start + input$numericMaxSearch)]
        # olap_dt[, ref_end := (ref_end - input$numericMaxSearch)]
        # olap_dt[, ref_width := (ref_width - 2*input$numericMaxSearch)]
        olap_dt[, ref_mid := (ref_end + ref_start)/2]
        
        olap_dt[, dist := min(abs(end - ref_mid), abs(start - ref_mid)), by = 1:nrow(olap_dt)]
        olap_dt[ref_mid < end & ref_mid > start, dist := 0]
        
        olap_dt[, o := rank(dist), by = ref_id]
        olap_dt = olap_dt[o <= input$sliderPerTSS,]
        rvOlapDt(olap_dt)
        enz_gr_sel = subset(enz_gr, name %in% olap_dt$name)
        rvEnzSelected(enz_gr_sel)
        
        tab = table(olap_dt[, .N, by = ref_id]$N)
        tab_num = as.numeric(tab)
        names(tab_num) = names(tab)
        p = ggplot() + geom_histogram(aes(x = names(tab_num), y = tab_num), stat = "identity") +
            labs(x = "Distance from promoter to nearest fragment", y = "count", title = input$selectCutter)
        
        p + cowplot::theme_cowplot()
    })
    
    output$plotOlapDist = renderPlot({
        req(rvAnnotExt())
        req(rvEnzSelected())
        annot_gr = rvAnnotExt()
        enz_gr = rvEnzSelected()
        dist_dt = as.data.table(distanceToNearest(annot_gr, enz_gr))
        ggplot(dist_dt[distance < 20000], aes(x = distance)) + 
            geom_histogram(bins = 50) + theme_cowplot()
    })
    
    output$ucscLink = renderUI({
        req(rvAnnotMissed())
        # req(rv$trackTxt)
        # tmpf = tempfile(pattern = "tracks_", tmpdir = CFG_DIR)
        # write.table(rv$trackTxt, file = tmpf, quote = F, row.names = F, col.names = F)
        tmp_url = sub(FILE_ROOT, URL_ROOT, rvConfigFile())
        
        ucsc_URL = "https://genome.ucsc.edu/cgi-bin/hgTracks?hgct_customText="
        ucsc_URL = paste0(ucsc_URL, tmp_url)
        ucsc_URL = paste0(ucsc_URL, "&guidelines=on/off")
        if(length(rvAnnotMissed()) == 0){
            pos = "chr21:5017493-5027492"
        }else{
            ###pos here
            rows <<- input$missedTssTable_state
            # message(paste(rows, collapse = ", "))
            save.image()
            pos = rvAnnotMissed()[1]
            strand(pos) = "*"
            pos = as.character(pos)    
        }
        pos = paste0("&position=", sub(":", "%3A", pos))
        ucsc_URL = paste0(ucsc_URL, pos)
        ucsc_URL = paste0(ucsc_URL, "&hgt.reset=1")
        ucsc_URL = paste0(ucsc_URL, "&db=", sub(",.+", "", input$selectGenome))
        tags$a(href = ucsc_URL, "to UCSC", target="_blank")
    })
    
    output$dlFrags = downloadHandler("frags.bed", 
                                     content = function(file){
                                         rtracklayer::export.bed(rvEnzSelected(), file)
                                     }
    )
    
    output$dlFragsPlusSeq = downloadHandler("fragsPlusSeq.bed", 
                                            content = function(file){
                                                library(BSgenome.Hsapiens.UCSC.hg38)
                                                library(Biostrings)
                                                library(data.table)
                                                
                                                gr = rvEnzSelected()
                                                seq = getSeq(Hsapiens, names = gr)
                                                gr$seq = seq
                                                df = data.frame(gr)
                                                # df$width = NULL
                                                df = df[, c(1:3, 6:7, 5, 8:9)]
                                                write.table(df, file = file, sep = "\t", row.names = F, col.names = F, quote = F)
                                                # rtracklayer::export.bed15(gr, format = "bed", file, expNames = c("CG", "seq"))
                                                # rtracklayer::export(gr, format = "bed", file, extraCols = c("CG" = "numeric", "seq" = "character"))
                                            }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

