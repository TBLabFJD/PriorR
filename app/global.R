library(shiny)
library(DT)    #  for data tables
library(dplyr)
library(shinyWidgets)
library(lazyeval)
library(shinydashboard)
library(shinydashboardPlus)
library(data.table)
library(shinyalert)
library(filesstrings)
library(shinyBS)
library(shinyjs)
library(optparse)


option_list=list(
  make_option(c('-p','--programdir'),type="character", help="Program directory."),
  make_option(c('-s','--sampledir'),type="character", help="Sample directory where you keep the diff files."))

opt <- parse_args(OptionParser(option_list = option_list))

dependencies_route <- paste(opt$programdir, 'Dependencies', sep='/')
Paneles_Orphanet <- paste(opt$programdir, 'Paneles_Orphanet', sep='/')
sample_route <- paste(opt$sampledir, 'DIFF', sep='/')

options(shiny.maxRequestSize = 1500*1024^2)

callback <- function(rows){
  c(
    sprintf("var rows = [%s];", toString(rows)),
    "$('#LOH').on('click', function(){",
    "    for(var i=0; i<rows.length; ++i){",
    "      var row = table.row(rows[i]);",
    "      if(row.length){",
    "        row.node().style.backgroundColor = ",
    "         $(this).prop('checked') ? 'yellow' : '';",
    "      }",
    "    }",
    "})"
  )
}

hpo_table <- paste(dependencies_route, "ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt",  sep='/')

hpos <- read.table(hpo_table, fill = TRUE, quote = "", 
                   sep = '\t', na.strings=c("",".","NA"), colClasses = NA)

# DOWNLOAD MODULE

downloadObjUI <- function(id) {
  ns <- NS(id)
  downloadButton(ns("data_download"), label = "Download Filtered Data", class = "btn-primary")}

downloadObj <- function(input, output, session, data) {
  output$data_download <- downloadHandler(
    filename = function() {
      paste("filetered_data", Sys.Date(), ".csv", sep = "")},
    content = function(file) {
      write.csv(data(), file)}
  )}