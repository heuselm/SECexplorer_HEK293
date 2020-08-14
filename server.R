# SECexplorer HEK293
# Exploring SEC-SWATH-MS data of the cycling HEK293 cells
# Heusel & Bludau et al. 2019 https://www.embopress.org/doi/10.15252/msb.20188438
#################################################################################
# Server functions
##################

## Load packages
library(shiny)
library(ggplot2)
library(CCprofiler)
library(plotly)

## Prepare data

# load data
protTraces <- readRDS("data/protTraces.rds")
calibration_functions = readRDS("data/calibration.rds")

# assign annotation columns to vectors
annotations <- names(protTraces$trace_annotation)
for (i in seq_along(annotations)){
  assign(annotations[i], protTraces$trace_annotation[[i]])
}

# Define server fundtions
shinyServer(function(input, output) {
  
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(get(input$fcolumn)))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Choose proteins/genes (type to search)",
                   values,
                   multiple = TRUE, 
                   options = list(maxOptions = 60000),
                   selected = sort(protTraces$trace_annotation$Gene_names[grep("COPS", protTraces$trace_annotation$Gene_names)]))
  })
  
  ## generate selected protein SEC traces plot
  output$plot <- renderPlotly({
    
    # Collect indices of the proteins matching the selection 
    indices <- numeric()
    for (i in 1:length(input$fvalue)){
      column_index <- which(annotations == input$fcolumn)
      indices = c(indices, grep(input$fvalue[i], get(input$fcolumn)))
    }
    
    # subset the data for these, including Mass column
    protTraces.s <- subset(protTraces, trace_subset_ids = protTraces$traces$id[indices])
    
    # Plot
    traces.long <- toLongFormat(protTraces.s$traces)
    traces.long <- merge(traces.long, protTraces.s$trace_annotation, by = "id", all.y = FALSE, all.x = TRUE)
    p <- ggplot(traces.long) +
      xlab('fraction') + ylab('Top2 peptide intensity sum') + 
      theme_bw() +
      theme(legend.position=input$legend.position) 
    if (is.null(calibration_functions)){
      q <- p + geom_line(aes(x=fraction, y=intensity, color=Gene_names))
      ggplotly(q)
    } else {
      lx.frc <- seq(10,(ncol(protTraces.s$traces)-1),10)
      lx <- paste( lx.frc , round(calibration_functions$FractionToMW(lx.frc), 1) , sep = '\n' )
      
      q <- p + geom_line(aes(x=fraction, y=intensity, color=Gene_names)) +
        scale_x_continuous(breaks = lx.frc , labels = lx) + xlab("fraction\n(apparent MW [kDa])")
      
      if(input$show_monomer_markers == TRUE){
        q <- q + geom_point(aes(x = calibration_functions$MWtoFraction(protein_mw), y = -10, color = Gene_names), pch = 23, fill = "white", size =2)
      }
      
      ggplotly(q)
    }
  })
  
  output$table <- DT::renderDataTable({
    # Collect indices of the proteins matching the selection 
    indices <- numeric()
    for (i in 1:length(input$fvalue)){
      column_index <- which(annotations == input$fcolumn)
      indices = c(indices, grep(input$fvalue[i], get(input$fcolumn)))
    }
    
    protTraces.s <- subset(protTraces, trace_subset_ids = protTraces$traces$id[indices])
    # and subset the data for these, including Mass column
    protTraces.s$trace_annotation[, id:=NULL]
  })
})
