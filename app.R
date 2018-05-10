#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
DEV=FALSE
library(shiny)
library(shinycssloaders)
library(BiocFileCache)
library(GenomicRanges)
library(data.table)
library(ggplot2)
ref_urls = list("hg38" = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz")
ref_files = list()
ref_bfc = BiocFileCache(cache = ".cache")
for(i in seq_along(ref_urls)){
    fpath = bfcrpath(x = ref_bfc, ref_urls[[i]], rnames = names(ref_urls)[i]) 
    ref_files[[names(ref_urls)[i]]] = fpath
}
if(DEV & !exists("tmp_ref")){
    tmp_ref = rtracklayer::import.gff(ref_files[[1]], format = "gtf", feature.type = "transcript")
    tmp_ref = subset(tmp_ref, seqnames == "chr21")
} 
enz_files = c("HindIII" = "~/HiC-Pro/data/frags_hg38canon_HindIII.bed",
              "DpnII" = "~/HiC-Pro/data/frags_hg38canon_MboI.bed")
if(DEV & !exists("tmp_enz")){
    tmp_enz = rtracklayer::import.bed(enz_files[1])
    tmp_enz = subset(tmp_enz, seqnames == "chr21")
} 
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("HiC capture probe design"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("selectGenome", "Genome", choices = "hg38"),
            tags$hr(),
            selectInput("selectCutter", "Enzyme", choices = c("DpnII", "HindIII"), selected = "HindIII"),
            withSpinner(plotOutput("plotFragSize")),
            tags$hr(),
            h5("TSS Extension"),
            fluidRow(
                column(numericInput("numericTSSdownstream", "Downstream", value = 5000, min = 0, max = 100000, step = 500), width = 6),
                column(numericInput("numericTSSupstream", "Upstream", value = 5000, min = 0, max = 100000, step = 500), width = 6)
            ),
            tags$hr(),
            numericInput("numericMaxSearch", "Max Search Extension", value = 10000, min = 0, max = 100000, step = 500),
            tags$hr(),
            h5("Allowed Fragment Sizes"),
            sliderInput("sliderFragSize", "Range", min = 0, max = 10000, step = 10, dragRange = TRUE, value = c(200, 6000)),
            sliderInput("sliderPerTSS", "Fragments per TSS", value = 1, min = 1, max = 20, step = 1),
            sliderInput("sliderPerFragment", "Probes per fragment", value = 2, min = 1, max = 10, step = 1),
            withSpinner(plotOutput("plotTSShits"))
            
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            actionButton("btnRestart", "Start Over"),
            actionButton("btnRemoveSelected", "Remove Selected"),
            actionButton("btnLimitToSelected", "Limit To Selected"),
            actionButton("btnRemoveVisible", "Remove Filtered"),
            actionButton("btnLimitToVisible", "Limit To Filtered"),
            withSpinner(DT::dataTableOutput(outputId = "filterTable")),
            tags$hr(),
            withSpinner(htmlOutput("calcMb")),
            tags$hr(),
            downloadButton("dlFrags", label = "Download Fragments")
            
        )
    )
)

server <- function(input, output) {
    
    rvTSScount = reactiveVal()
    rvGeneAnnotRaw = reactiveVal()
    rvGeneAnnotFiltering = reactiveVal()
    rvGeneAnnotExt = reactiveVal()
    # rvGeneAnnotFinal = reactiveVal()
    
    rvEnz = reactiveVal()
    
    rvEnzSelected = reactiveVal()
    
    #when gene ref is finalized and enz is loaded,
    #find nearest n enz to each item in gene ref
    # observeEvent(
    #     eventExpr = {
    #         rvGeneAnnotExt()
    #         rvEnz()
    #     }, 
    #     handlerExpr = {
    #         ref_gr = rvGeneAnnotExt()
    #         mcols(ref_gr) = NULL
    #         ref_gr = unique(ref_gr)
    #         ref_gr$id = paste0("tss_", seq_along(ref_gr))
    #         end(ref_gr) = end(ref_gr) + input$numericMaxSearch
    #         start(ref_gr) = start(ref_gr) - input$numericMaxSearch
    #         # ref_gr = resize(ref_gr, 1, "center")
    #         enz_gr = rvEnz()
    #         
    #         olaps = as.data.table(findOverlaps(ref_gr, enz_gr, ignore.strand = TRUE))
    #         # fine = olaps[, .N, by = queryHits][N <= input$sliderPerTSS]$queryHits
    #         
    #         olap_dt = cbind(as.data.table(ref_gr[olaps$queryHits]), as.data.table(enz_gr[olaps$subjectHits]))
    #         nc = ncol(as.data.frame(head(ref_gr)))
    #         colnames(olap_dt)[1:nc] = paste0("ref_", colnames(olap_dt)[1:nc])
    #         
    #         olap_dt[, ref_start := (ref_start + input$numericMaxSearch)]
    #         olap_dt[, ref_end := (ref_end - input$numericMaxSearch)]
    #         olap_dt[, ref_mid := (ref_end + ref_start)/2]
    #         
    #         olap_dt[, dist := min(abs(end - ref_mid), abs(start - ref_mid)), by = 1:nrow(olap_dt)]
    #         olap_dt[ref_mid < end & ref_mid > start, dist := 0]
    #         olap_dt[, o := order(dist, decreasing = FALSE), by = ref_id]
    #         olap_dt = olap_dt[o <= input$sliderPerTSS,]
    #         
    #         # tab = table(olap_dt[, .N, by = ref_id]$N)
    #         # tab_num = as.numeric(tab)
    #         # names(tab_num) = names(tab)
    #         # p = ggplot() + geom_histogram(aes(x = names(tab_num), y = tab_num), stat = "identity") +
    #         #     labs(x = "number of fragments per tss", y = "count", title = input$selectCutter)
    #         
    #         enz_gr_sel = subset(enz_gr, name %in% olap_dt$name)
    #         rvEnzSelected(enz_gr_sel)
    #         # hist(olaps[, .N, by = queryHits]$N)
    #     }
    # )
    
    # observeEvent(
    #     eventExpr = {
    #         input$selectCutter
    #     }, 
    #     handlerExpr = {
    #         # print("go")
    #         fpath = enz_files[input$selectCutter]
    #         gr = rtracklayer::import.bed(fpath)
    #         rvEnz(gr)
    #     })
    # 
    #if selected reference changes, update rvGeneAnnotRaw
    observeEvent(
        eventExpr = {
            input$selectGenome
        }, 
        handlerExpr = {
            fpath = ref_files[[input$selectGenome]]
            if(DEV){
                rvGeneAnnotRaw(as.data.frame(tmp_ref))
            }else{
                rvGeneAnnotRaw(as.data.frame(rtracklayer::import.gff(fpath, format = "gtf", feature.type = "transcript")))
            }
            
            
        }
    )
    
    #if rvGeneAnnotRaw changes or start over, init rvGeneAnnotFiltering
    observeEvent(
        eventExpr = {
            input$btnRestart
            rvGeneAnnotRaw()
        }, 
        handlerExpr = {
            rvGeneAnnotFiltering(rvGeneAnnotRaw())
        }
    )
    
    observeEvent(input$btnRemoveSelected, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        
        rsel = sort(as.integer(input$filterTable_rows_selected))
        rvis = setdiff(rvis, rsel)
        if(length(rsel) < 1){
            showNotification("Empty selection.", type = "error")
        }else{
            rvGeneAnnotFiltering(rvGeneAnnotFiltering()[rvis,])
            showNotification(paste("Removed", length(rsel), "entries."))
        }
    })
    
    observeEvent(input$btnLimitToSelected, {
        rsel = sort(as.integer(input$filterTable_rows_selected))
        nr = nrow(rvGeneAnnotFiltering())
        if(length(rsel) < 1){
            showNotification("Empty selection.", type = "error")
        }else{
            rvGeneAnnotFiltering(rvGeneAnnotFiltering()[rsel,])
            showNotification(paste("Removed", nr - length(rsel), "entries."))
        }
    })
    
    observeEvent(input$btnRemoveVisible, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        nr = nrow(rvGeneAnnotFiltering())
        if(length(rvis) < 1){
            showNotification("Would be empty.", type = "error")
        }else if(length(rvis) == nr){
            showNotification("No filters.", type = "error")
        }else{
            rvGeneAnnotFiltering(rvGeneAnnotFiltering()[-rvis,])
            showNotification(paste("Removed", length(rvis), "entries."))
        }
        
    })
    
    observeEvent(input$btnLimitToVisible, {
        rvis = sort(as.integer(input$filterTable_rows_all))
        nr = nrow(rvGeneAnnotFiltering())
        if(length(rvis) < 1){
            showNotification("Would be empty.", type = "error")
        }else if(length(rvis) == nr){
            showNotification("No filters.", type = "error")
        }else{
            rvGeneAnnotFiltering(rvGeneAnnotFiltering()[rvis,])
            showNotification(paste("Removed", nr - length(rvis), "entries."))
        }
        
    })
    
    #if TSS size parameter or rvGeneAnnotFiltering change, update rvGeneAnnotExt
    observeEvent(
        eventExpr = 
        {
            rvGeneAnnotFiltering()
            input$numericTSSdownstream
            input$numericTSSupstream
        }, 
        handlerExpr = {
            gr = GRanges(rvGeneAnnotFiltering())
            gr = promoters(gr, 
                           upstream = input$numericTSSupstream, 
                           downstream = input$numericTSSdownstream)
            
            rvGeneAnnotExt(gr)
        }
    )
    
    # #as rvGeneAnnotFiltering or merge paramters change, update rvGeneAnnotFinal
    # observeEvent({
    #     
    # })
    
    output$filterTable = DT::renderDataTable({
        DT::datatable(rvGeneAnnotFiltering(), filter = "top", options = list(pageLength = 10, scrollX = T)) #%>%
    }) 
    
    output$calcMb = renderUI({
        req(rvGeneAnnotFiltering())
        req(rvEnzSelected())
        rvis = sort(as.integer(input$filterTable_rows_all))
        ntss = nrow(rvGeneAnnotFiltering()[rvis,])
        total = 120 * input$sliderPerFragment * input$sliderPerTSS * ntss
        total_str = round(total / 10^6, digits = 2)
        
        nr_enz = length(rvEnzSelected())
        nr_total = nr_enz*120 * input$sliderPerFragment
        nr_total_str = round(nr_total / 10^6, digits = 2)
        tags$p(
            tags$span(paste("TSS Count:", ntss)),
            tags$br(),
            tags$span(paste("Fragments Per TSS:", input$sliderPerTSS)),
            tags$br(),
            tags$span(paste("Probes Per Fragment:", input$sliderPerTSS)),
            tags$br(),
            tags$span("Probe Size: 120 bp"),
            tags$br(),
            tags$span(paste("Possible Total:", total_str, "Mb")),
            tags$hr(),
            tags$span(paste("non-redundant Fragments:", nr_enz)),
            tags$br(),
            tags$span(paste("non-redundant Total:", nr_total_str, "Mb"))
            
        )
    })
    
    output$plotFragSize = renderPlot({
        req(input$selectCutter)
        input$selectCutter
        fpath = enz_files[input$selectCutter]
        if(DEV){
            gr = tmp_enz    
        }else{
            gr = rtracklayer::import.bed(fpath)    
        }
        
        
        rvEnz(gr)
        
        df = data.frame(width = width(gr))
        qmax = quantile(df$width, .999)
        p = ggplot(df[df$width < qmax,, drop = FALSE], aes(x = width)) + geom_histogram(bins = 50)
        bp = ggplot2::ggplot_build(p)
        xrng = bp$layout$panel_ranges[[1]]$x.range
        yrng = bp$layout$panel_ranges[[1]]$y.range
        p + annotate("text", x = mean(xrng), y = mean(yrng), hjust = 0, label = paste0("mean ", round(mean(df$width)), 
                                                                                       "\nmedian ", median(df$width))) +
            labs(title = "Fragment Size Distribution") + cowplot::theme_cowplot()
    })
    
    output$plotTSShits = renderPlot({
        req(rvGeneAnnotExt())
        req(rvEnz())
        ref_gr = rvGeneAnnotExt()
        mcols(ref_gr) = NULL
        ref_gr = unique(ref_gr)
        ref_gr$id = paste0("tss_", seq_along(ref_gr))
        end(ref_gr) = end(ref_gr) + input$numericMaxSearch
        start(ref_gr) = start(ref_gr) - input$numericMaxSearch
        enz_gr = rvEnz()
        
        frag_rng = input$sliderFragSize
        k = width(enz_gr) >= min(frag_rng) & width(enz_gr) <= max(frag_rng)
        enz_gr = enz_gr[k]
        
        olaps = as.data.table(findOverlaps(ref_gr, enz_gr, ignore.strand = TRUE))
        
        olap_dt = cbind(as.data.table(ref_gr[olaps$queryHits]), as.data.table(enz_gr[olaps$subjectHits]))
        nc = ncol(as.data.frame(head(ref_gr)))
        colnames(olap_dt)[1:nc] = paste0("ref_", colnames(olap_dt)[1:nc])
        
        olap_dt[, ref_start := (ref_start + input$numericMaxSearch)]
        olap_dt[, ref_end := (ref_end - input$numericMaxSearch)]
        olap_dt[, ref_mid := (ref_end + ref_start)/2]
        
        olap_dt[, dist := min(abs(end - ref_mid), abs(start - ref_mid)), by = 1:nrow(olap_dt)]
        olap_dt[ref_mid < end & ref_mid > start, dist := 0]
        olap_dt[, o := order(dist, decreasing = FALSE), by = ref_id]
        olap_dt = olap_dt[o <= input$sliderPerTSS,]
        
        enz_gr_sel = subset(enz_gr, name %in% olap_dt$name)
        rvEnzSelected(enz_gr_sel)
        
        tab = table(olap_dt[, .N, by = ref_id]$N)
        tab_num = as.numeric(tab)
        names(tab_num) = names(tab)
        p = ggplot() + geom_histogram(aes(x = names(tab_num), y = tab_num), stat = "identity") +
            labs(x = "number of fragments per tss", y = "count", title = input$selectCutter)
        
        p + cowplot::theme_cowplot()
    })
    
    output$dlFrags = downloadHandler("frags.bed", 
                                     content = function(file){
                                         rtracklayer::export.bed(rvEnzSelected(), file)
                                     }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

