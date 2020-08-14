# SECexplorer HEK293
# Exploring SEC-SWATH-MS data of the cycling HEK293 cells
# Heusel & Bludau et al. 2019 https://www.embopress.org/doi/10.15252/msb.20188438
#################################################################################
# User Interface
################

## Load packages
library(shiny)
library(ggplot2)
library(CCprofiler)
library(plotly)
library(DT)

## Prepare data

# load data
protTraces <- readRDS("data/protTraces.rds")
calibration_functions = readRDS("data/calibration.rds")

# assign annotation columns to vectors
annotations <- names(protTraces$trace_annotation)
for (i in seq_along(annotations)){
  assign(annotations[i], protTraces$trace_annotation[[i]])
}

## Define User Interface 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("SECexplorer HEK293"),
  
  # Sidebar to select proteins for display
  sidebarPanel(
    selectInput("fcolumn", "Choose annotation type for protein/gene selection", annotations, selected = "Gene_names"),
    
    uiOutput("fcolumnvalues"), 
    #The possible choices of this field are calculated on the server side and passed over by uiOutput
    
    checkboxInput("show_monomer_markers", label = "Show monomer expected elution fraction markers", value = TRUE),
    
    textInput("legend.position", label = "Plot legend position (or \"none\")", value = "right"),
    
    p("Select genes/protein of interest above to display elution profiles."),
    p("Reference:\n \nHeusel & Bludau et. al. Mol Syst Biol. 2019 Jan 14;15(1):e8438. doi: 10.15252/msb.20188438"),
    p("Contact:\naebersold@imsb.biol.ethz.ch\nmoritz.heusel@med.lu.se")
  ),
  
  # mainPanel: Show interactive line plot the selected Traces
  mainPanel(
    plotlyOutput("plot"),
    DT::dataTableOutput("table")
  )
))
