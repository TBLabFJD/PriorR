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
library(shinyFiles)
library(plotly)
library(ggplot2)


programdir='/session/app'  

dependencies_route <- paste(programdir, 'Dependencies', sep='/')
annotation_route <- paste(programdir, 'resources', sep='/')
www_route <- paste(programdir, 'www', sep='/')
logo_route <- paste(www_route, 'logo.png', sep='/')
Paneles_Orphanet <- paste0(programdir, '/','Paneles_Orphanet', '/')  

options(shiny.maxRequestSize = 2500*1024^2)

callback <- function(rows, rowsr, rowso){
  c(
    sprintf("var rows = [%s];", toString(rows)),
    "$('#ROH').on('click', function(){",
    "    for(var i=0; i<rows.length; ++i){",
    "      var row = table.row(rows[i]);",
    "      if(row.length){",
    "        row.node().style.backgroundColor = ",
    "         $(this).prop('checked') ? 'salmon' : '';",
    "      }",
    "    }",
    "})",
    sprintf("var rowsr = [%s];", toString(rowsr)),
    "$('#prior').on('click', function(){",
    "    for(var i=0; i<rowsr.length; ++i){",
    "      var row = table.row(rowsr[i]);",
    "      if(row.length){",
    "        row.node().style.backgroundColor = ",
    "         $(this).prop('checked') ? 'red' : '';",
    "      }",
    "    }",
    "})",
    sprintf("var rowso = [%s];", toString(rowso)),
    "$('#prior').on('click', function(){",
    "    for(var i=0; i<rowso.length; ++i){",
    "      var row = table.row(rowso[i]);",
    "      if(row.length){",
    "        row.node().style.backgroundColor = ",
    "         $(this).prop('checked') ? 'darkorange' : '';",
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


### DASHBOARD DESIGN ###

ui <- function(request){
  
  dashboardPagePlus( skin="green",
                     
                     ##############
                     ### HEADER ###
                     ##############
                     
                     header= dashboardHeaderPlus(  
                       title = "PriorR 3.0",
                       left_menu = tagList (
                         dropdownBlock(
                           id="Help",
                           title="",
                           icon="question-circle",
                           actionButton("pdf", "Help", onclick = "window.open('manual.pdf')", style = "width:100px"))
                       )),
                     
                     ###############
                     ### SIDEBAR ###
                     ###############
                     
                     dashboardSidebar(disable = FALSE,
                                      width =250,
                                      
                                      ########## SNVs ##########
                                      
                                      # INPUT #
                                      
                                      div(id='tab1_sidebar',
                                          
                                          # SIDEBAR #
                                          
                                          sidebarMenu(id = "sidebarmenu",
                                                      
                                                      # NAME #
                                                      useShinyalert(),  # Set up shinyalert
                                                      verbatimTextOutput("name", placeholder = FALSE),
                                                      actionButton("preview", "  Sign in", icon=icon("user"), width = '200px', style="color: #fff; background-color: #9bcd9b; border-color: #f0f0f0 ; font-size: 17px"),
                                                      # UPLOAD DATA #
                                                      menuItem("Upload analysis data",
                                                               tabName = "upload",
                                                               icon=icon("upload"),
                                                               fileInput("file1", "Upload your SNV File",
                                                                         multiple = FALSE,
                                                                         accept = c("text/csv",
                                                                                    "text/comma-separated-values,text/plain", ".tsv")),
                                                               actionBttn(
                                                                 inputId = "WGS",
                                                                 label = "Analysis of WGS", 
                                                                 style = "minimal",
                                                                 color = "success"
                                                               )),
                                                      
                                                      
                                                      # COLUMNS #
                                                      menuItem("Columns",
                                                               tabName = "columns",
                                                               icon=icon("columns"),
                                                               pickerInput("Columns",   
                                                                           choices = NULL,
                                                                           selected = NULL,
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE )),
                                                      # FILTERS #
                                                      
                                                      menuItem("Filters",
                                                               tabName = "Filters",
                                                               icon=icon("ruler-horizontal"),
                                                               
                                                               pickerInput("Type", "Type:",   
                                                                           choices = NULL,
                                                                           selected =  NULL,
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE ),
                                                               
                                                               pickerInput("Region", "Region:",   
                                                                           choices = c("3UTR" ,"5UTR" , "DOWNSTREAM", "EXONIC" ,"INTRONIC" ,"ncRNA","regulatory","SPLICING","UPSTREAM", NA ),
                                                                           selected =  c("EXONIC", "SPLICING"),
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE ),
                                                               
                                                               
                                                               pickerInput("Consequence","Consequence:",
                                                                           choices = NULL,
                                                                           selected = NULL,
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE ),
                                                               
                                                               
                                                               fileInput("file2", "Regions",
                                                                         multiple = FALSE,
                                                                         accept = c("text/csv",
                                                                                    "text/comma-separated-values,text/plain",
                                                                                    ".csv")),
                                                               prettyCheckbox(inputId= "CANONICAL", label = "CANONICAL", value = TRUE, 
                                                                              outline= TRUE,  bigger = TRUE, status = 'success', width = NULL),
                                                               
                                                               prettyCheckbox(inputId="predicted", label = "EXCLUDE PREDICTED", value = FALSE, 
                                                                              outline= TRUE,  status = 'success', width = NULL)),
                                                      
                                                      # FREQUENCIES #
                                                      
                                                      menuItem(
                                                        "Frequecies",
                                                        setSliderColor(c("ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen"), c(1, 2, 3, 4, 5, 6, 7, 8)),
                                                        icon=icon("users"),
                                                        sliderInput("freq1", "gnomADg AF", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq2", "gnomADg AF popmax", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq3", "gnomADe AF", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq5", "FJD_MAF AF", min = 0, max = 1, value = 0.01)),
                                                      
                                                      
                                                      # PATHOGENICITY # 
                                                      
                                                      menuItem(
                                                        "Pathogenicity",
                                                        setSliderColor("ForestGreen", 1),
                                                        icon=icon("file-medical-alt"),
                                                        sliderTextInput("patho1", "CADD PHRED", choices = seq(from = 0, to = 75, by = 5), selected = 15), 
                                                        prettyCheckbox(inputId= "Clinvar", label = "NO Clinvar's Benign Variants", value = T,  outline= TRUE,  bigger = TRUE, status = 'success', width = NULL)),
                                                      
                                                      
                                                      # INHERITANCE #
                                                      
                                                      menuItem(
                                                        "Inheritance",           
                                                        tabName = "Inheritance",
                                                        icon=icon("child"),
                                                        pickerInput("Inheritance", "Inheritance Pattern:",
                                                                    choices= c("","Dominant",
                                                                               "Recessive-Homozygous",
                                                                               "Recessive-Compound heterozygous",
                                                                               "X-linked"), selected = "", multiple = FALSE ),
                                                        fileInput("file11", "For a family-based study upoload a pedegree file",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv"))
                                                      ),
                                                      
                                                      # GENES #
                                                      
                                                      menuItem(
                                                        "genes",
                                                        tabName = "genes",
                                                        icon=icon("dna"),
                                                        textInput("text", "Your genes"),
                                                        fileInput("file3", "Your Virtual Panel",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv", ".txt")),
                                                        actionBttn('reset1', 'Reset', style = "bordered", color='success'),
                                                        fileInput("file4", "Exclude Virtual Panel",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv", ".txt")),
                                                        actionBttn('reset2', 'Reset', style = "bordered", color='success'),
                                                        selectizeInput("orphanet", "Orphanet panel & Glow Genes prioritization",
                                                                       multiple=TRUE,
                                                                       selected = NULL, 
                                                                       options =list(placeholder='Please select a panel'),
                                                                       before_last_dot(list.files(Paneles_Orphanet, pattern='_panel_glow.txt$'))),                                              
                                                        sliderInput("glow", "GLOW genes", min = -1, max = 1800, value = -1)
                                                        
                                                      ),
                                                      
                                                      # HPO #
                                                      
                                                      menuItem(
                                                        "Phenotypes / HPO",
                                                        tabName = "Phenotype",
                                                        icon=icon("user-cog"),
                                                        selectizeInput("Phenotype",
                                                                       "Phenotype",
                                                                       multiple=TRUE,
                                                                       selected = NULL, 
                                                                       options = list(placeholder='Please select a phenotype'),
                                                                       choices = levels(hpos$V3))),
                                                      
                                                      # ROH #
                                                      
                                                      menuItem(
                                                        "Family info / ROH",
                                                        tabName = "ROH",
                                                        icon=icon("home"),
                                                        textInput("text2", "Probandus", width = '100%'),
                                                        fileInput("file5", "Pedegree File",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv")),
                                                        prettyCheckbox(inputId="ROH", label = "ROH", value = FALSE, outline= TRUE,  status = 'success', width = NULL)
                                                      ),
                                                      
                                                      
                                                      
                                                      # COMMERCIAL PIPELINE  #
                                                      
                                                      menuItem(
                                                        "Commercial Pipeline",       
                                                        tabName = "Commercial Pipeline",
                                                        icon=icon("project-diagram"),
                                                        fileInput("file10", "Upload commercial vcf",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv", ".vcf")),
                                                        
                                                        prettyCheckbox(inputId= "not-commercial", label = "COLOURING", value = F, 
                                                                       outline= TRUE,  bigger = TRUE, status = 'success', width = NULL),
                                                        
                                                        prettyCheckbox(inputId="not-commercial-filtering", label = "FILTERING", value = F, 
                                                                       outline= TRUE,  status = 'success', width = NULL)
                                                        
                                                      ),
                                                      
                                                      
                                                      # CLASSIFICATION AND ANNOTATION #
                                                      
                                                      menuItem(
                                                        "Classification", 
                                                        icon=icon("th-list"),
                                                        prettyCheckbox(inputId="Classification", label = "Annotate Classified Variants", value = FALSE,outline= TRUE, fill = TRUE, status = 'success', width = NULL)),
                                                      
                                                      menuItem(
                                                        "Prioritization", 
                                                        icon=icon("list-alt"),
                                                        # actionBttn(inputId = 'prior', 'Prioritization', style = "bordered", color='success')),
                                                        
                                                        prettyCheckbox(inputId="prior", label = "COLORING", value = FALSE, outline= TRUE, fill = TRUE, status = 'success', width = NULL),
                                                        # actionBttn('prior1', 'Prioritization', style = "bordered", color='success')),
                                                        
                                                        prettyCheckbox(inputId="prior1", label = "PRIORITIZATION", value = FALSE, outline= TRUE, fill = TRUE, status = 'success', width = NULL)),
                                                      
                                                      
                                                      menuItem(
                                                        "Save Data",
                                                        icon=icon("save"),
                                                        downloadObjUI(id = "download1"))
                                                      
                                          )),
                                      
                                      ########## CNVs #########
                                      
                                      shinyjs::hidden(
                                        
                                        div(id='tab2_sidebar',
                                            fileInput("file6", "Upload your CNV File",
                                                      multiple = FALSE,
                                                      accept = c("text/csv",
                                                                 "text/comma-separated-values,text/plain",
                                                                 ".csv", ".tsv")),
                                            textInput("text5", "Sample"),
                                            sidebarMenu(id = "sidebarmenu2",
                                                        
                                                        # COLUMNS #
                                                        
                                                        menuItem("Columns",
                                                                 tabName = "Columns_cnv",
                                                                 icon=icon("columns"),
                                                                 pickerInput("Columns_cnv",   
                                                                             choices = NULL,
                                                                             selected = NULL,
                                                                             options = list(`actions-box` = TRUE),
                                                                             multiple = TRUE )),
                                                        
                                                        
                                                        # OPTIONS #
                                                        menuItem("Filters",
                                                                 setSliderColor(c("ForestGreen"), 1),
                                                                 tabName = "Filters",
                                                                 icon=icon("ruler-horizontal"),
                                                                 
                                                                 pickerInput("Annotation_mode", "Annotation Mode:",   
                                                                             choices = c("full", "split"),
                                                                             selected =  "full",
                                                                             options = list(`actions-box` = TRUE),
                                                                             multiple = TRUE ),
                                                                 
                                                                 
                                                                 pickerInput("Type_cnv", "Type:",   
                                                                             choices = NULL,
                                                                             selected =  NULL,
                                                                             options = list(`actions-box` = TRUE),
                                                                             multiple = TRUE ),
                                                                 
                                                                 pickerInput("Region_cnv", "Region:",   
                                                                             choices = NULL,
                                                                             selected =  NULL,
                                                                             options = list(`actions-box` = TRUE),
                                                                             multiple = TRUE ),
                                                                 
                                                                 
                                                                 radioGroupButtons("Programs","Number of programs",
                                                                                   choices = c('1', '2', '3'),
                                                                                   selected = '3',
                                                                                   status = "success"),
                                                                 
                                                                 sliderTextInput(
                                                                   inputId = "acmg",
                                                                   label = "ACMG class:", 
                                                                   choices = seq(from=5, to=1, by=-1),
                                                                   selected = 3
                                                                 )
                                                                 
                                                                 
                                                        ),
                                                        
                                                        
                                                        # GENE ANNOTATION #
                                                        
                                                        menuItem(
                                                          "Genes",
                                                          tabName = "genes",
                                                          icon=icon("dna"),
                                                          textInput("text6", "Your genes"),
                                                          fileInput("file8", "Your Virtual Panel",
                                                                    multiple = FALSE,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain",
                                                                               ".csv")),
                                                          fileInput("file9", "Exclude Virtual Panel",
                                                                    multiple = FALSE,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain",
                                                                               ".csv")),
                                                          sliderInput("glow_cnv", "GLOW genes", min = 0, max = 1800, value = 0),
                                                          selectizeInput("orphanet_cnv", "Gene Panel Orphanet",
                                                                         multiple=TRUE,
                                                                         selected = NULL, 
                                                                         options =list(placeholder='Please select a panel'),
                                                                         before_last_dot(list.files(Paneles_Orphanet, pattern='_panel.txt$')))),
                                                        
                                                        # CLASSIFICATION #
                                                        
                                                        menuItem(
                                                          "Classification",
                                                          icon=icon("th-list"),
                                                          prettyCheckbox(inputId="Classification", label = "Annotate Classified Variants", value = FALSE, outline= TRUE, fill = TRUE, status = 'success', width = NULL)),
                                                        
                                                        menuItem(
                                                          "Frequencies",
                                                          setSliderColor(c("ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen"), c(7, 8, 9, 10)),
                                                          icon=icon("users"),
                                                          sliderInput("freq1_cnv", "Gain AF max", min = 0, max = 1, value = 0.01),
                                                          sliderInput("freq2_cnv", "Loss AF max", min = 0, max = 1, value = 0.01)),
                                                        
                                                        menuItem(
                                                          "Save Data",
                                                          icon=icon("save"),
                                                          downloadObjUI(id = "download2"))
                                                        
                                                        # 
                                                        
                                            )
                                        ))
                     ),
                     
                     ############
                     ### BODY ###
                     ############
                     
                     dashboardBody(
                       useShinyjs(),
                       #                    tags$head(tags$style(HTML('
                       #   .modal.in .modal-dialog{
                       #     width:100%;
                       #     height:100%;
                       #     margin:0px;
                       #   }
                       # 
                       #   .modal-content{
                       #     width:100%;
                       #     height:100%;
                       #   }
                       # '))),
                       tabsetPanel(
                         id = "page",
                         type = "hidden",
                         
                         #Your landing page
                         tabPanelBody("landing-page",
                                      tags$style("#land2 {font-size:48px;
                                        color:black;
                                        display:block; }"),
                                      div(
                                        style = "position: absolute;
                   left: 0;
                   top: 0;
                   z-index: 10000;
                   width: 100%;
                   height: 100%;
                   background: lightblue;",
                                        div(
                                          p(
                                            fluidRow(
                                              column(1,imageOutput("land1")),
                                            ),
                                            fluidRow(
                                              column(5, offset = 4, htmlOutput('land2')),
                                            ), 
                                            fluidRow(
                                              column(5, offset = 3, actionButton("close-landing-page", "Take me to Priorr"),
                                                     column(5, offset = 3, actionBttn(
                                                       inputId = "annotate",
                                                       label = "Annotate my VCF",
                                                       style = "minimal")),
                                              )
                                            )
                                          ),
                                          #                      style = "position: relative;
                                          # top: 40%;
                                          # left: 40%;",
                                          # h1("Welcome to Priorr"),
                                          # #Button to close landing page
                                          # actionBttn(
                                          #   inputId = "annotate",
                                          #   label = "Annotate my VCF", 
                                          #   style = "minimal"
                                          # ),
                                          # actionButton("close-landing-page", "Take me to Priorr"),
                                          bsModal("survey1", "Select annotation information","annotate", 
                                                  
                                                  awesomeRadio(inputId = "assembly", label = "Choose assembly", choices = c("GRCh37", "GRCh38"), selected = NULL, inline = TRUE, checkbox = TRUE),
                                                  fileInput("file14", "Upload your VCF to annotate", multiple = FALSE, accept = c(".vcf", ".vcf.gz")),
                                                  actionButton("EnterVCF", "Annotate file"),
                                                  uiOutput("downloadvcf"))
                                                  
                                        )
                                      )
                                      
                         ),
                         tabPanelBody('body-content',
                                      tabsetPanel(
                                        id="navbar",
                                        tabPanel(title="SNV",id="tab1",value='tab1_val', DTOutput("contents")),
                                        tabPanel(title="CNV",id="tab2",value='tab2_val', DTOutput("contents_cnv"))
                                      ),
                                      tags$head(tags$style(HTML('
        .main-header .logo {
        background-color: #ad1d28;
        font-family:"Georgia";
        font-weight: bold;
        font-size: 52 px;
        }
        .navbar{
        background-color: cyan
        }
        '
                                      ))),
                                      tags$style("#info1 {font-size:20px;
                                       color:black;
                                       display:block; }"),
                                      # tags$style("#info2 {font-size:14px;
                                      #             display:block;
                                      #            }"),
                                      
                                      tags$style("#info4 {font-size:16px;
                                       color:green;
                                       display:block; }"),
                                      tags$style("#info5 {font-size:13px;
                                       display:block;
                                      }"),
                                      tags$style("#info6 {font-size:16px;
                                       color:green;
                                       display:block; }"),
                                      tags$style("#info7 {font-size:13px;
                                       display:block;
                                      }"),
                                      tags$style("#info8 {font-size:16px;
                                       color:green;
                                       display:block; }"),
                                      tags$style("#info9 {font-size:13px;
                                       display:block;
                                      }"),
                                      div(id='tab1_body',
                                          
                                          boxPlus(
                                            title = "Variant Information", 
                                            closable = FALSE, 
                                            width = NULL,
                                            status = "success", 
                                            solidHeader = FALSE, 
                                            collapsible = TRUE,
                                            collapsed = TRUE,
                                            
                                            p(
                                              fluidRow(
                                                column(8, textOutput('info1')),
                                                column(3, htmlOutput('info4'))
                                              ), 
                                              fluidRow(
                                                column(4, verbatimTextOutput('info2')),
                                                column(3, offset = 0.5,
                                                       HTML("<div style='height: 20px;'>"),
                                                       plotOutput("info3", width = 230, height = 150),
                                                       HTML("</div>"),
                                                ),
                                                column(3,  offset=1, htmlOutput('info5'))
                                              ),
                                              fluidRow(
                                                column(8, textOutput('info6')),
                                                column(3, textOutput('info8')),
                                              ),
                                              fluidRow(
                                                column(4, htmlOutput('info7')),
                                                column(3, offset = 0.5,
                                                       HTML("<div style='height: 20px;'>"),
                                                       plotOutput("info10", width = 235, height = 125),
                                                       HTML("</div>"),
                                                ),
                                                column(3, offset=1, htmlOutput('info9')),
                                              )
                                            ))),
                                      boxPlus(
                                        title = "Saved Sessions", 
                                        closable = FALSE, 
                                        width = NULL,
                                        status = "warning", 
                                        solidHeader = FALSE, 
                                        collapsible = TRUE,
                                        collapsed = TRUE,
                                        
                                        p(
                                          fluidRow(
                                            column(3, textInput(inputId = "description", label = "Bookmark description")), 
                                            column(4, bookmarkButton(id="bookmarkBtn"))
                                          ),
                                          dataTableOutput("urlTable", width = "100%"))
                                      ),
                                      
                                      # SURVEY WGS #
                                      
                                      bsModal("survey", "Select WGS data information","WGS", 
                                              prettyCheckbox(inputId="canonical_filters", label = "Canonical", value = TRUE, outline= TRUE, fill = TRUE, status = 'success', width = NULL),
                                              awesomeRadio("freq_filters", "Frequency filters (Gnomadg AF)", choices=list("rare (0.01)"=0.01, "very rare (0.001)"=0.001, "ultra-rare (0.0001)"=0.0001), selected = NULL, inline=T, status = 'success'),
                                              
                                              fileInput("file12", "Panel", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")),
                                              prettyCheckbox(inputId="codding", label = "Codding regions", outline= TRUE, fill = TRUE, status = 'success', width = NULL, value = TRUE,  bigger = TRUE),
                                              fileInput("file13", "WGS variant file", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".tsv", ".txt"))),
                                      
                                      # SURVEY ANNOTATION #
                                      
                                      
                                      
                                      dataTableOutput("hpos_table"),
                                      uiOutput("modal"),
                                      tags$script("Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
Shiny.onInputChange(variableName, null);
});
")
                         )
                       )            
                     ))
  
}

# Define server logic to read selected file ----
server <- function(input, output, session) {
  
  
  ##################
  ### ANNOTATION ###
  ##################
  
  observeEvent(ignoreInit=T, input$file14,{ 
    
    rv$vcfdownload <- NULL
    
    file = input$file14
    ext = tools::file_ext(file$datapath)
    
    sharedir = annotation_route
    assembly = input$assembly
    inputvcf = input$file14$datapath
    
    withProgress(message = 'Annotating...', value = 1/2, {
      validate(need(ext %in% c("vcf", "vcf.gz"), "Please upload a VCF file"))
      sample_prefix = gsub("(.vcf|.vcf.gz)$", "", file$name, perl = T)
      rv$vcfdownload <- paste0(sample_prefix, ".pvm.tsv")
      system(paste('bash /session/app/annotation.sh', sharedir, assembly, inputvcf))
        })
    
    })
  
  
  ####################
  ### UPLOAD FILES ###
  ####################
  
  rv <- reactiveValues(df=NULL)
  
  # VARIANTS
  
  df1 <- reactive({1
    if(!is.null(input$file1)){
      df <- read.table(input$file1$datapath, fill = TRUE, quote = "", header = TRUE,
                       sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
    }
    
  })
  
  
  df2 <- reactive({
    if(!is.null(input$file13)){
      
      if (input$canonical_filters ==T){canonical='-c'}
      if (input$codding == T){codding ='-r'}
      if (!is.null(input$file12)) {panel= paste('-p', input$file12$datapath)} else {panel=NULL}
      if (length(input$file13$datapath)!=0){
        withProgress(message='Filtering WGS data...',   {
          system(paste('/opt/conda/envs/Priorr/bin/python3.10 /session/app/scripts/wgs_filters.py', '-f', input$file13$datapath, canonical, '-af', input$freq_filters, codding, panel) )} ) }
      
      
      df <- read.table('/session/app/tmp/wgs.tsv', fill = TRUE, quote = "", header = TRUE,
                       sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
    }
  })
  
  
  
  
  
  ###  condition this observer to display df1()
  observeEvent(df1(), {
    rv$df <- df1()
  })
  
  ###  condition this observer to display df2()
  observeEvent(df2(), {
    rv$df <- df2()
  })
  
  
  df_cnv <- reactive({
    req(input$file6)
    df_cnv <- read.delim(input$file6$datapath, fill = TRUE, header = TRUE,
                         sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
    df_cnv <- df_cnv %>% mutate( ACMG_class= gsub("full=",  "", as.character(ACMG_class) ))
    print(head(df_cnv))
    # df_cnv$ACMG_class = gsub("full=",  "", as.character(df_cnv$ACMG_class) )
  })
  
  
  
  
  # PANEL
  
  panel <- reactive({
    req(input$file3)
    rv1$data <- read.csv(input$file3$datapath,
                         sep = ',', header=FALSE)
  })
  
  excluded_panel <- reactive({
    req(input$file4)
    rv2$data <- read.csv(input$file4$datapath,
                         sep = ',', header=FALSE)
  })
  
  panel_cnv <- reactive({
    req(input$file8)
    panel <- read.csv(input$file8$datapath,
                      sep = ',', header=FALSE)
  })
  
  excluded_panel_cnv <- reactive({
    req(input$file9)
    panel <- read.csv(input$file9$datapath,
                      sep = ',', header=FALSE)
  })
  
  # BED FILE
  
  bed <- reactive({
    req(input$file2)
    bed <- read.csv(input$file2$datapath,
                    sep = '\t', header=FALSE)
  })
  
  # cnvs
  
  bed <- reactive({
    req(input$file7)
    bed_cnv <- read.csv(input$file7$datapath,
                        sep = '\t', header=FALSE)
  })
  
  # PED FILE
  
  ped <- reactive({
    req(input$file5)
    ped <- read.csv(input$file5$datapath,
                    sep = '\t', header=FALSE)
  })
  
  ped_inh <- reactive({
    req(input$file11)
    ped_inh <- read.csv(input$file11$datapath,
                        sep = '\t', header=FALSE)
  })
  
  
  # COMMERCIAL VCF
  
  comm <- reactive({
    req(input$file10)
    comm <- read.csv(input$file10$datapath,
                     sep = '\t', header=FALSE, comment.char = '#')
  })
  
  
  fjd_variants <- reactive({fjd_variants <- diff()[diff()["IN_FILE"] == 1,]})
  
  
  # RESET
  
  rv1 <- reactiveValues(
    data = NULL,
    clear = FALSE
  )
  
  rv2 <- reactiveValues(
    data = NULL,
    clear = FALSE
  )
  
  
  
  #################
  ### FUNCTIONS ###
  #################
  
  ##############
  # USER LOGIN #
  ##############
  observeEvent(input$preview, {
    shinyalert("User Name", type = "input") # Show a modal when the button is pressed
  })
  
  
  output$name <- renderText( paste("User:", input$shinyalert),
                             outputArgs = list(
                             ))
  
  
  observeEvent( input$shinyalert, {
    dir.create(file.path(programdir, input$shinyalert))
    setwd(file.path(programdir, input$shinyalert))
  })
  
  getwd()
  
  #############
  # FUNCTIONS #
  #############
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  
  #############################################
  ### FILTERS, CHOICES, FREQUENCIES, PANELS ###
  #############################################
  
  ################
  # FILTERS SNVs #
  ################
  
  
  observeEvent(rv$df, {
    req(rv$df)
    vf <- glob2rx("*VF")
    gt <- glob2rx("*GT")
    ad <- glob2rx("*AD")
    dp <- glob2rx("*DP")
    updatePickerInput(session, inputId = "Columns", choices = colnames(rv$df), 
                      selected = c('CHROM', 'POS', 'REF', 'ALT', 'Location', 'SYMBOL', 'Gen_Full_name', 'VARIANT_CLASS',  grep(gt, colnames(rv$df), value = TRUE), grep(vf, colnames(rv$df), value = TRUE), grep(ad, colnames(rv$df), value = TRUE), grep(dp, colnames(rv$df), value = TRUE), 'Existing_variation', 'Genomic_region', 'CANONICAL', 'Feature','BIOTYPE','Consequence', 'HGVSc', 'HGVSp', 'CLNSIG', 'OMIM_phenotype', 'Orphanet_disorder','CADD_PHRED', 'gnomADg_AF', 'gnomADg_AF_popmax', 'gnomADg_popmax', 'gnomADe_AF', 'FJD_MAF_AF', 'CSVS_AF', 'MutScore','ada_score', 'rf_score', 'N_Pathogenic_pred', 'N_Benign_pred', 'N_predictions', 'kaviar_AF', 'ExACpLI'))
    updatePickerInput(session, inputId = "Consequence", choices = levels(addNA(rv$df$Consequence)), selected = levels(addNA(rv$df$Consequence)) )
    updatePickerInput(session, inputId = "Type", choices = levels(as.factor(rv$df$VARIANT_CLASS)), selected = levels(as.factor(rv$df$VARIANT_CLASS)) )
    
    print(typeof(rv$df$VARIANT_CLASS))
  })
  
  
  observeEvent(input$file3, {
    rv1$clear <- FALSE
    print(rv1$clear)
  })
  
  observeEvent(input$file4, {
    rv2$clear <- FALSE
    print(rv2$clear)
  })
  
  filtered_df <- reactive( {
    
    DF <- rv$df %>%
      
      select(input$Columns) %>%
      filter(Consequence %in% input$Consequence )  %>%
      filter(Genomic_region %in% input$Region ) %>%
      filter(VARIANT_CLASS  %in% input$Type) %>%
      filter(if (input$CANONICAL == FALSE) input$CANONICAL==input$CANONICAL else CANONICAL == 'YES') %>%
      filter(if (input$predicted == FALSE) input$predicted == input$predicted  else !grepl("X", Feature)) %>%
      filter(if (input$text == '') SYMBOL == SYMBOL else SYMBOL %in% list(input$text) ) %>%
      filter( as.numeric(gnomADg_AF) <= input$freq1 | is.na(gnomADg_AF) ) %>%
      filter( as.numeric(gnomADg_AF_popmax) <= input$freq2 | is.na(gnomADg_AF_popmax) ) %>%
      filter( as.numeric(gnomADe_AF) <= input$freq3 | is.na(gnomADe_AF))  %>%
      filter( as.numeric(FJD_MAF_AF) <= input$freq5 | is.na(FJD_MAF_AF))  %>%
      filter( CADD_PHRED >= input$patho1 | is.na(CADD_PHRED))
    
    # Inheritance pattern
    
    gt <- glob2rx("*GT")
    GT_column <- grep(gt, colnames(rv$df), value = TRUE)
    
    print(GT_column)
    
    if(is.null(input$file11) ){
      if(input$Inheritance == 'Dominant'){
        DF <- DF[(DF[[GT_column]] == '0/1' | DF[[GT_column]] == '0|1' ) & (grepl("dominant", DF[,"OMIM_phenotype"], ignore.case=TRUE )),]
      }
      else if(input$Inheritance == 'Recessive-Homozygous'){
        DF <- DF[(DF[[GT_column]] == '1/1' | DF[[GT_column]] == '1|1' | DF[[GT_column]] == '1/2' ) & (grepl("recessive", DF[,"OMIM_phenotype"], ignore.case=TRUE )) ,]
      }
      else if(input$Inheritance == 'Recessive-Compound heterozygous'  ){
        DF <- DF %>% group_by(SYMBOL) %>% filter(n()>1)
        DF <- DF[(DF[[GT_column]] == '0/1' | DF[[GT_column]] == '0|1' ) & (grepl("recessive", DF[,"OMIM_phenotype"], ignore.case=TRUE )) ,]
        
      }
      else if(input$Inheritance == 'X-linked'){
        DF <- DF  %>% filter(CHROM == 'chrX')
      }
    }
    else{
      affected = ped_inh()[ped_inh()$V6 == 2, 2]
      samples1 = lapply(affected, gsub, pattern ="-", replacement=".")
      affected = paste("X", samples1, "_GT",sep="")
      GT_affected = unlist(affected)
      print(GT_affected)
      
      not_affected = ped_inh()[ped_inh()$V6 == 1, 2]
      samples2 = lapply(not_affected, gsub, pattern ="-", replacement=".")
      not_affected = paste("X", samples2, "_GT", sep="")
      GT_not_affected = unlist(not_affected)
      
      father = levels(ped_inh()$V3)[2]
      father = gsub(pattern ="-", replacement=".", father)
      father_gt = paste0("X", father, "_GT")
      mother = levels(ped_inh()$V4)[2]
      mother = gsub(pattern ="-", replacement=".", mother)
      mother_gt = paste0("X", mother, "_GT")
      
      rel_not_affected = GT_not_affected[!(GT_not_affected %in% c(mother, father))]
      
      males_affected = ped_inh()[(ped_inh()$V5 == 1) & (ped_inh()$V6 == 2), 2]
      males_affected = gsub(pattern ="-", replacement=".", males_affected)
      males_affected = paste0("X", males_affected, "_GT")
      
      
      if(input$Inheritance == 'Dominant'){
        
        print('hello I am dominant')
        
        DF <- DF[which((rowSums((DF[GT_affected] =='0/1') | (DF[GT_affected] =='0|1'))  == length(GT_affected))
                       & (rowSums((DF[GT_not_affected]=='0/0') | (DF[GT_not_affected]=='0|0')) == length(GT_not_affected)) & grepl("dominant", DF$OMIM_phenotype, ignore.case=TRUE )),]
      } 
      else if(input$Inheritance == 'Recessive-Homozygous'){
        
        DF <- DF[which(rowSums((DF[GT_affected] =='1/1') | (DF[GT_affected] =='1|1'))  == length(GT_affected) & ((DF[father_gt]=='0/1') | (DF[father_gt]=='0|1') | (is.na(father_gt))) & ((DF[mother_gt]=='0/1') | (DF[mother_gt]=='0|1')| (is.na(mother_gt))) & rowSums((DF[rel_not_affected]!='1/1') | (DF[rel_not_affected]=='1|1') ) == length(rel_not_affected) & grepl("recessive", DF$OMIM_phenotype, ignore.case=TRUE )),]
      }
      
      else if (input$Inheritance == 'Recessive-Heterozygous') {
        DF <- DF %>% group_by(SYMBOL) %>% filter(n()>1)
        # DF <- DF[which(rowSums((DF[GT_affected] =='0/1') | (DF[GT_affected] =='0|1'))  == length(GT_affected) & grepl("recessive", DF$OMIM_phenotype, ignore.case=TRUE )),]
      }
      else if(input$Inheritance == 'X-linked'){
        DF <- DF  %>% filter(CHROM == 'chrX')
        DF <- DF[which(rowSums((DF[GT_affected] !='0/0') & (DF[GT_affected] =='0|0'))  == length(GT_affected) & rowSums((DF[males] =='1/1') | (DF[males_affected] =='1|1'))  == length(males_affected)),]
      }
    }
    
    if( input$glow != -1 ){
      DF <- DF  %>% filter(as.numeric(GLOWgenes) <= input$glow)
    }

    # Gene panels
    
    if (!is.null(input$file3) && rv1$clear == F) {
      DF <- DF %>% filter(SYMBOL %in%panel()[,1] )
      print(paste0('Value of clear before resetting is: ', rv1$clear))
    } 
    else if(!is.null(input$file3) && rv1$clear == T){
      DF$SYMBOL = DF$SYMBOL
    }
    
    # Exclude panel
    
    if(!is.null(input$file4) && rv2$clear == F){
      excluded_genes = unlist(excluded_panel()[,1])
      DF <- DF %>% filter(!SYMBOL %in% excluded_genes )
      print(excluded_genes)
      # DF <- DF %>% filter(!SYMBOL %in% excluded_panel()[,1] )
    }
    else if(!is.null(input$file4) && rv2$clear == T){
      DF$SYMBOL = DF$SYMBOL
    }
    
    # GLOW genes
    
    if (!is.null(input$orphanet)){
      print(input$orphanet)
      orphanet_path <- paste0(Paneles_Orphanet, input$orphanet, "_panel_glow.txt")
      print(orphanet_path)
      glow_file <- read.table( orphanet_path, header=FALSE, sep='\t' )
      colnames(glow_file) = c("SYMBOL", "V2", "GLOWgenes")
      DF = left_join(df, glow_file[, c("SYMBOL", 'GLOWgenes')], by='SYMBOL')
      DF = arrange(DF, GLOWgenes)
      
      # DF <- DF %>% filter(SYMBOL %in% gene_set[,1] )
    }
    
    # Bed file upload
    if (!is.null(input$file2)){
      DF$POS <- as.numeric(as.character(DF$POS))
      DF <- left_join( DF, bed(), c("CHROM" = "V1" )) %>% filter( POS >= V2, POS <= V3 )
    }
    
    # Classified variants annotation
    if (input$Classification == TRUE){
      variants_db <- read.table('/session/app/Database/variants_priorr.txt', fill = TRUE, 
                                header = TRUE,sep = '\t', na.strings=c("",".","NA"), colClasses = NA )
      variants_db$POS <- as.integer(variants_db$POS)
      # print(variants_db)
      DF <- left_join( DF, variants_db, by = c("CHROM", "POS", "REF", "ALT"))
    }
    
    # PHENOTYPE
    
    if (!is.null(input$Phenotype)){
      set <- hpos %>% filter( hpos$V3 %in% input$Phenotype)
      DF <- DF %>% filter( SYMBOL %in% set$V2 )
    }
    
    # Differencial variants
    
    if(!is.null(input$file10) && input$'not-commercial-filtering' == TRUE){
      DF <- DF[!(DF$CHROM %in% comm()$V1) & !(DF$POS %in% comm()$V2), ]
    }
    
    # Variant Classification
    DF$Classification.Button = shinyInput( actionButton, nrow(DF), 'button_', label = "Classify" , onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    
    if (input$Clinvar == TRUE) {
      DF <- DF %>% filter(!grepl("benign", CLNSIG, ignore.case=TRUE ))
    }
    
    # Prioritization
    
    if (input$prior1 == TRUE){
      print('here I am')
      # #   
      # #   #####
      # #   
      orange_class = which(
        (
          (as.numeric(DF[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(DF[,'gnomADg_AF_popmax'])>0.001)
          & (((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60) | (as.numeric(DF[,'N_predictions'])==0))
          & (as.numeric(DF[,'CADD_PHRED'])>25)
          & (as.numeric(DF[,'MutScore'])>0.80 | is.na(DF[,'MutScore']))
        ) |
          (
            (as.numeric(DF[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(DF[,'gnomADg_AF_popmax'])>0.001)
            & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60)
            & ((as.numeric(DF[,'CADD_PHRED'])>25) | (as.numeric(DF[,'MutScore'])>0.80))
          )  |
          (
            ( grepl('splice', DF[,'Consequence']) | grepl('intron', DF[,'Consequence']))
            & (as.numeric(DF[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(DF[,'gnomADg_AF_popmax'])>0.001)
            & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60)
            & ((DF[,'ada_score']>0.8) | (DF[,'rf_score']>0.8))
          ) |
          (
            (grepl('pathogenic', DF[, 'CLNSIG'],ignore.case=TRUE))
            & (DF[, 'CLNSIG'] != 'Conflicting_interpretations_of_pathogenicity')
            & (DF[, 'CLNSIG'] != 'Pathogenic')
            & (DF[, 'CLNSIG'] != 'Pathogenic/Likely_pathogenic')
            & (DF[, 'CLNSIG'] != 'Likely_pathogenic')
          ) |
          (
            (as.numeric(DF[,'gnomADg_AF_popmax'])<0.001 | is.na(DF[,'gnomADg_AF_popmax']))
            & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60)
            & ((as.numeric(DF[,'CADD_PHRED'])>25) | (as.numeric(DF[,'MutScore'])>0.80))
            & (DF[, 'CLNSIG'] == 'Uncertain_significance')
          )
        
      ) 
      print(orange_class)
      patho_conseq = c('start_lost', 'stop_gained', 'frameshift_variant', 'stop_lost', 'inframe_deletion', 'inframe_insertion')
      
      red_class = which(
        (
          (DF[, 'CLNSIG'] == 'Pathogenic') | (DF[, 'CLNSIG'] == 'Pathogenic/Likely_pathogenic') |(DF[, 'CLNSIG'] == 'Likely_pathogenic') |
            (
              (as.numeric(DF[,'gnomADg_AF_popmax'])<0.001 | is.na(DF[,'gnomADg_AF_popmax']))
              & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60)
              & (as.numeric(DF[,'CADD_PHRED'])>25)
              & (as.numeric(DF[,'MutScore'])>0.80 | is.na(DF[,'MutScore']))
            ) |
            (
              (as.numeric(DF[,'gnomADg_AF_popmax'])<0.001 | is.na(DF[,'gnomADg_AF_popmax']))
              & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.60)
              & (as.numeric(DF[,'CADD_PHRED'])>25 | is.na(DF[,'CADD_PHRED']))
              & (as.numeric(DF[,'MutScore'])>0.80)
            ) |
            (
              (DF[,'Consequence'] %in% patho_conseq)
              & (as.numeric(DF[,'gnomADg_AF_popmax'])<0.001 | is.na(DF[,'gnomADg_AF_popmax']))
              & ((as.numeric(DF[,'N_Pathogenic_pred'])/as.numeric(DF[,'N_predictions']))>0.75)
              & (is.na(DF[,'CADD_PHRED']) & is.na(DF[,'MutScore']))
            ) |
            (
              (grepl('splice', DF[,'Consequence']) | grepl('intron', DF[,'Consequence']))
              & (as.numeric(DF[,'gnomADg_AF_popmax'])<0.001 | is.na(DF[,'gnomADg_AF_popmax']))
              & (as.numeric(DF[,'ada_score'])>0.8 | as.numeric(DF[,'rf_score'])>0.8)
            )
        ) &
          !(
            grepl("benign", DF[, 'CLNSIG'], ignore.case=TRUE )
          ) &
          !(
            grepl("Uncertain_significance", DF[, 'CLNSIG'], ignore.case=TRUE )
          )
      )
      
      print(red_class)
      #   #####
      # #   
      rows <- 1:nrow(DF)
      rest <- rows[!(rows %in% c(red_class, orange_class))]
      prior_order <- c(red_class, orange_class, rest)
      print(prior_order)
      DF <- DF[prior_order,]
    }
    # # 
    return(DF)
  })
  
  
  observeEvent(input$reset1, {
    rv1$data <- NULL
    rv1$clear <- TRUE
    reset('file3')
  })
  
  observeEvent(input$reset2, {
    rv2$data <- NULL
    rv2$clear <- TRUE
    reset('file4')
  })
  
  
  ################
  # FILTERS CNVs #
  ################
  
  
  observeEvent(df_cnv(), {
    req(df_cnv())
    updatePickerInput(session, inputId = "Columns_cnv", choices = colnames(df_cnv()), selected = c( "ACMG_class","Samples_ID", "SV_chrom", "SV_start", "SV_end",  "SV_length","SV_type", "N_PROGRAMS", "Annotation_mode", "Gene_name", "Gene_count", "Tx" ,"Frameshift", "Exon_count", "Location", "Location2", "P_gain_phen" , "P_gain_hpo", "P_loss_phen",  "P_loss_hpo",  "P_snvindel_phen", "B_gain_AFmax", "B_loss_AFmax", "DDD_disease", "OMIM_phenotype"))
    updatePickerInput(session, inputId = "Region_cnv", choices = levels(df_cnv()$Location2), selected = levels(df_cnv()$Location2)) 
    updatePickerInput(session, inputId = "Type_cnv", choices = levels(as.factor(df_cnv()$SV_type)), selected = levels(as.factor(df_cnv()$SV_type))) 
  })
  
  
  filtered_df_cnv <- reactive( {
    
    DF <- df_cnv() %>%
      
      select(input$Columns_cnv) %>%
      filter(if (input$text5 == '') Samples_ID == Samples_ID else Samples_ID %in% list(input$text5) ) %>%
      filter(Annotation_mode %in% input$Annotation_mode ) %>%
      filter(SV_type %in% input$Type_cnv ) %>%
      filter(ACMG_class >= input$acmg ) %>%
      filter(N_PROGRAMS %in% input$Programs ) %>%
      filter(if (input$text6 == '') Gene_name == Gene_name else grepl(input$text6, Gene_name)) %>%
      filter(if (is.null(input$file8)) Gene_name == Gene_name else grepl(paste( panel_cnv()[,1], collapse= '|'), Gene_name )) %>%
      filter(if (is.null(input$file9)) Gene_name == Gene_name else !grepl(paste( excluded_panel_cnv()[,1], collapse= '|'), Gene_name )) %>% 
      filter( as.numeric(B_gain_AFmax) <= input$freq1_cnv | is.na(B_gain_AFmax))  %>%
      filter( as.numeric(B_loss_AFmax) <= input$freq2_cnv | is.na(B_loss_AFmax)) 
    
    
    # # GLOW genes
    # 
    # if( input$glow_cnv != 0 ){
    #   DF <- DF  %>% filter(as.numeric(GLOWgenes) <= input$glow_cnv)
    # }
    # 
    
    
    # Orphanet gene sets
    # 
    # if (!is.null(input$orphanet_cnv)){
    #   orphanet_path <- paste0(Paneles_Orphanet, input$orphanet_cnv, ".txt")
    #   gene_set <- read.table( orphanet_path, header=FALSE, sep='\t' )
    #   # print(gene_set)
    #   DF <- DF %>% filter(grepl(paste( gene_set[,1], collapse= '|'), Gene_name ))
    # }
    # #
    # 
    # 
    
    
    
    
    # Orphanet gene sets
    
    # Bed file upload
    # if (!is.null(input$file7)){
    #   DF$START <- as.numeric(as.character(DF$POS))
    #   DF <- left_join( DF, bed_cnv(), c("CHROM" = "V1" )) %>% filter( POS >= V2, POS <= V3 )
    # }
    # 
    # Classified variants annotation
    # if (input$Classification == TRUE){
    #   variants_db <- read.table('/home/rrr/n.txt', fill = TRUE, 
    #                             header = TRUE,sep = '\t', na.strings=c("",".","NA"), colClasses = NA )
    #   variants_db$POS <- as.character(variants_db$POS)
    #   # print(variants_db)
    #   DF <- left_join( df(), variants_db, by = c("CHROM", "POS", "REF", "ALT"))
    # }
    # 
    # 
    # # Variant Classification
    # DF$Classification.Button = shinyInput( actionButton, nrow(DF), 'button_', label = "Classify" , onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    # 
    
  })
  
  #######################
  ### OTHER FUNCTIONS ###
  #######################
  
  ####################################
  # MODAL CLASSIFICATION POP_UP SNVs #
  ####################################
  
  observeEvent(input$select_button, {
    s <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    output$modal <- renderUI({
      tagList(
        bsModal(paste('model', s ,sep=''), "Classification", "select_button", size = "small",
                radioButtons("pathogenicity", "pathogenicity", choices=c("Benign", "Possibly benign", "Unknown", "Pathogenic","Possibly pathogenic"), selected = NULL),
                textAreaInput("comments", label = h3("Comments") , value = "", width = "100%", height = "100px", resize = "none"),
                actionButton("Enter", "Save")
        ))
    })
    toggleModal(session,paste('model', s ,sep=''), toggle = "Assessment")
    
  })      
  
  ####################################
  # ANNOTATION OF VARIANTS TO SNV DB #
  ####################################
  
  observeEvent(input$Enter, {
    selectedRow = as.numeric(strsplit(input$select_button, "_")[[1]][2])
    line = paste( filtered_df()[selectedRow,1],  filtered_df()[selectedRow,2],  filtered_df()[selectedRow,3], filtered_df()[selectedRow,4], input$pathogenicity, input$comments, input$shinyalert, sep='\t')
    print(line)
    write(line, '/session/app/Database/variants_priorr.txt', append = TRUE)
    
  })
  
  
  #######
  # ROH #
  #######
  
  yellowRows <- reactive({
    
    req(filtered_df())
    
    
    # Only probandus
    if( input$text2 != ""){
      samples = input$text2
      ROH = paste( samples, "ROH", sep= "_")
      print("probandus is:")
      print(ROH)
      ind = names(filtered_df()) %in% ROH
      which(filtered_df()[, ind] == "True")  - 1L
      
      # With ped file
    }else if (!is.null(input$file5)){
      affected = ped()[ped()$V6 == 2, 2]
      samples1 = lapply(affected, gsub, pattern ="-", replacement=".")
      samples = paste("X", samples1, sep="")
      
      print(paste0("Affected individuals are: ",affected))
      
      
      if(length(samples) > 1){
        print("Here I am")
        ROH = lapply(samples, paste, "ROH", sep = "_")
        ind = names(filtered_df()) %in% ROH
        which(rowSums(filtered_df()[, ind] == "True") == sum(ind))  - 1L
        
      }else{
        ROH = paste( samples, "ROH", sep= "_")
        ind = names(filtered_df()) %in% ROH
        which(filtered_df()[, ind] == "True")  - 1L}
      
      # Single sample  (only one _ROH column)
    }else {
      sample_ROH = names(filtered_df())[endsWith(names(filtered_df()), "_ROH")]
      # print("No ped file needed, the only sample is:")
      # print(sample_ROH)
      if(length(sample_ROH) == 1 ){
        samples = strsplit(sample_ROH, "_")[[1]][1]
        ROH = paste( samples, "ROH", sep= "_")
        ind = names(filtered_df()) %in% ROH
        which(filtered_df()[, ind] == "True")  - 1L
      }else{
        return(NULL)}
    }
    
  })
  
  print(yellowRows)
  
  ##########################
  # VARIANT PRIORITIZATION #
  ##########################
  
  red_class <- reactive({
    req(filtered_df())
    patho_conseq = c('start_lost', 'stop_gained', 'frameshift_variant', 'stop_lost', 'inframe_deletion', 'inframe_insertion')
    
    red_class = which(
      (
        (filtered_df()[, 'CLNSIG'] == 'Pathogenic') | (filtered_df()[, 'CLNSIG'] == 'Pathogenic/Likely_pathogenic') |(filtered_df()[, 'CLNSIG'] == 'Likely_pathogenic') |
          (
            (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.001 | is.na(filtered_df()[,'gnomADg_AF_popmax'])) 
            & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60)
            & (as.numeric(filtered_df()[,'CADD_PHRED'])>25) 
            & (as.numeric(filtered_df()[,'MutScore'])>0.80 | is.na(filtered_df()[,'MutScore'])) 
          ) | 
          (
            (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.001 | is.na(filtered_df()[,'gnomADg_AF_popmax']))
            & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60)
            & (as.numeric(filtered_df()[,'CADD_PHRED'])>25 | is.na(filtered_df()[,'CADD_PHRED']))
            & (as.numeric(filtered_df()[,'MutScore'])>0.80)
          ) |
          (
            (filtered_df()[,'Consequence'] %in% patho_conseq)
            & (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.001 | is.na(filtered_df()[,'gnomADg_AF_popmax']))
            & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.75)
            & (is.na(filtered_df()[,'CADD_PHRED']) & is.na(filtered_df()[,'MutScore']))
          ) |
          (
            (grepl('splice', filtered_df()[,'Consequence']) | grepl('intron', filtered_df()[,'Consequence']))
            & (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.001 | is.na(filtered_df()[,'gnomADg_AF_popmax']))
            & (as.numeric(filtered_df()[,'ada_score'])>0.8 | as.numeric(filtered_df()[,'rf_score'])>0.8)
          )
      ) &
        !(
          grepl("benign", filtered_df()[, 'CLNSIG'], ignore.case=TRUE )
        ) &       
        !(
          grepl("Uncertain_significance", filtered_df()[, 'CLNSIG'], ignore.case=TRUE )
        ) 
    )- 1L
    
    print(red_class)
    
  })
  
  
  orange_class <- reactive({
    req(filtered_df())
    
    orange_class = which(
      (
        (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(filtered_df()[,'gnomADg_AF_popmax'])>0.001)
        & (((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60) | (as.numeric(filtered_df()[,'N_predictions'])==0))
        & (as.numeric(filtered_df()[,'CADD_PHRED'])>25)
        & (as.numeric(filtered_df()[,'MutScore'])>0.80 | is.na(filtered_df()[,'MutScore'])) 
      ) |
        (
          (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(filtered_df()[,'gnomADg_AF_popmax'])>0.001)
          & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60) 
          & ((as.numeric(filtered_df()[,'CADD_PHRED'])>25) | (as.numeric(filtered_df()[,'MutScore'])>0.80))
        )  |
        (
          ( grepl('splice', filtered_df()[,'Consequence']) | grepl('intron', filtered_df()[,'Consequence'])) 
          & (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.01 &  as.numeric(filtered_df()[,'gnomADg_AF_popmax'])>0.001) 
          & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60)  
          & ((filtered_df()[,'ada_score']>0.8) | (filtered_df()[,'rf_score']>0.8))
        ) |
        (
          (grepl('pathogenic', filtered_df()[, 'CLNSIG'],ignore.case=TRUE)) 
          & (filtered_df()[, 'CLNSIG'] != 'Conflicting_interpretations_of_pathogenicity')
          & (filtered_df()[, 'CLNSIG'] != 'Pathogenic') 
          & (filtered_df()[, 'CLNSIG'] != 'Pathogenic/Likely_pathogenic') 
          & (filtered_df()[, 'CLNSIG'] != 'Likely_pathogenic')
        ) |
        (
          (as.numeric(filtered_df()[,'gnomADg_AF_popmax'])<0.001 | is.na(filtered_df()[,'gnomADg_AF_popmax'])) 
          & ((as.numeric(filtered_df()[,'N_Pathogenic_pred'])/as.numeric(filtered_df()[,'N_predictions']))>0.60)
          & ((as.numeric(filtered_df()[,'CADD_PHRED'])>25) | (as.numeric(filtered_df()[,'MutScore'])>0.80))
          & (filtered_df()[, 'CLNSIG'] == 'Uncertain_significance') 
        )   
      
    ) - 1L
    
    print(orange_class)
  })
  
  
  
  ###################
  ### TOGGLE MENU ###
  ###################
  
  
  values <- reactiveValues(selectedTab = 1)
  
  observeEvent( input$navbar, {
    toggle("tab1_sidebar", condition = input$navbar == "tab1_val")
    toggle("tab2_sidebar", condition = input$navbar == "tab2_val")
    toggle("tab1_body", condition = input$navbar == "tab1_val")
  })
  
  
  ### OUTPUT TABLE ###
  
  FJD_rows <- reactive({
    if (!is.null(input$file10) && input$'not-commercial' == TRUE){
      FJD_rows = !(filtered_df()$CHROM %in% comm()$V1) & !(filtered_df()$POS %in% comm()$V2)
    }else if (is.null(input$file10)){
      FJD = FALSE
    }
    print(FJD_rows)
    return(FJD_rows)
  })
  
  ### OUTPUT TABLE ###
  
  output$contents <- renderDT({
    req(filtered_df())
    datatable(
      filtered_df(),
      filter = "top", 
      class = "display nowrap compact",
      editable=T,
      rownames = F,
      selection = 'single',
      callback = JS(callback(yellowRows(),red_class(), orange_class())), 
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = '40px', targets ="_all" )),
        pageLength = 10),
      escape = FALSE)},
    server = FALSE)
  
  
  output$contents_cnv <- renderDT({
    req(filtered_df_cnv())
    datatable(
      filtered_df_cnv(),
      class = "display nowrap compact", 
      filter = "top", 
      selection = 'single',
      options = list(
        scrollX = TRUE, columnDefs = list(list(autoWidth = TRUE, width = '10px', targets = c(3))),
        pageLength = 10),      
      escape = FALSE) %>%
      formatStyle( "ACMG_class",
                   backgroundColor = styleEqual( c(1,2,3,4,5), c("chartreuse3", "greenyellow", "yellow", "orange", "red" )))
  },
  
  server = FALSE)
  
  
  ####################
  ### LANDING PAGE ###
  ####################
  
  
  output$land1 = renderImage({
    filename <- normalizePath(file.path(logo_route))
    list(src = filename,
         contentType = 'image/png',
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  
  output$land2 = renderPrint({
    cat('Welcome to PriorR')
  })
  
  ####################
  ### VARIANT INFO ###
  ####################
  
  output$info1 = renderPrint({
    s = input$contents_rows_selected
    if (length(s)) {
      cat('VARIANT')
    }
  })
  
  output$info2 = renderPrint({
    s = input$contents_rows_selected
    GT_column <- grep('_GT', colnames(rv$df), value = TRUE)
    samples = gsub('_GT', '', GT_column)
    variant = paste(c(as.character(filtered_df()[s, 'CHROM']), ' ', as.character(filtered_df()[s, 'POS']), ' ', as.character(filtered_df()[s, 'REF']), ' ', as.character(filtered_df()[s, 'ALT'])), collapse = '-')
    gt_info=c()
    for (gt in GT_column){
      sample = gsub('_GT', '', gt)
      gt_info =c(gt_info, sample, as.character(filtered_df()[s, gt]))
    }
    gt_info = paste(gt_info, collapse= '  ')
    clinvar_info = paste("CLINVAR:", as.character(filtered_df()[s, 'CLNSIG']), sep='')
    # variant_info = paste( variant, as.character(filtered_df()[s, GT_column]), collapse = '\t')
    # variant_vector = c( variant_info,  as.character(filtered_df()[s, 'SYMBOL']), as.character(filtered_df()[s, 'Consequence']), as.character(filtered_df()[s, 'HGVSc']), as.character(filtered_df()[s, 'HGVSp']))
    if (length(s)) {
      cat(variant, '\n')
      cat(gt_info, '\n')
      cat(as.character(filtered_df()[s, 'SYMBOL']), '\n')
      cat(as.character(filtered_df()[s, 'Consequence']), '\n')
      cat(as.character(filtered_df()[s, 'HGVSc']), '\n')
      cat(as.character(filtered_df()[s, 'HGVSp']), '\n')
      cat(clinvar_info, '\n')
    }
  })
  
  output$info3 = renderPlot({ 
    
    print(colors)
    s = input$contents_rows_selected
    dp <- glob2rx("*DP")
    DP_column <- grep(dp, colnames(rv$df), value = TRUE)
    palette = c("#1B9E77","#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
    colors = palette[1:length(DP_column)]
    print(colors)
    
    df_filtered <- rv$df %>% distinct(CHROM, POS, .keep_all = T)
    DP_variant <-  as.numeric(filtered_df()[s, DP_column])
    
    sample_column =c()
    dp_sample_column =c()
    for (s in DP_column){
      sample =  gsub('_DP', '', s)      
      sample_rep = rep(sample, nrow(df_filtered))
      sample_column =c(sample_column, sample_rep)
      dp_sample_column = c(dp_sample_column, df_filtered[, s])
    }
    
    DP_df = data.frame(sample=sample_column, DP= dp_sample_column)
    
    DP_df[, 'DP' ] = sapply(DP_df[, 'DP' ], as.numeric)
    ggplot(DP_df, aes(DP, color=sample)) +
      geom_density() +
      scale_x_continuous(trans='log10') +
      theme_classic()  + 
      geom_vline(xintercept = DP_variant, color = colors) +
      scale_color_brewer(palette="Dark2")
    
  }, height = 125, width = 300)
  
  output$info4 = renderPrint({
    s = input$contents_rows_selected
    if (length(s)) {
      cat('Population Frequency')
    }
  })
  
  output$info5 = renderUI({
    s = input$contents_rows_selected
    if (length(s)) {
      pop_vector = c('FJD MAF:', as.character(filtered_df()[s, 'FJD_MAF_AF']), 'Gnomadg MAF:', as.character(filtered_df()[s, 'gnomADg_AF']), 'Spanish Frequency:', as.character(filtered_df()[s, 'CSVS_AF']))
      HTML(paste(pop_vector, collapse = "<br/>"))
      # HTML(paste('FJD MAF:', as.character(filtered_df()[s, 'FJD_MAF_AF']), sep="<br/>"))
      # HTML(paste('Gnomadg MAF:', as.character(filtered_df()[s, 'gnomADg_AF']), sep="<br/>"))
      # HTML(paste('Spanish Frequency:', as.character(filtered_df()[s, 'CSVS_AF']), sep="<br/>"))
    }
  })
  
  output$info6 = renderPrint({
    s = input$contents_rows_selected
    if (length(s)) {
      cat('OMIM Phenotype:')
    }
  })
  
  output$info7 = renderUI({
    s = input$contents_rows_selected
    if (length(s)) {
      # phen_list = strsplit(as.character(filtered_df()[s, 'OMIM_phenotype']), ";" )
      phen_list = gsub(";", "<br/n>", as.character(filtered_df()[s, 'OMIM_phenotype']))
      # print(phen_list)
      HTML(paste(phen_list, sep="<br/>"))
    }
  })
  
  output$info8 = renderPrint({
    s = input$contents_rows_selected
    if (length(s)) {
      cat('Pathogenicity predictors:')
    }
  })
  
  output$info9 = renderUI({
    s = input$contents_rows_selected
    if (length(s)) {
      patho_vector = c('CADD Phred:', as.character(filtered_df()[s, 'CADD_PHRED']), 'Mut-Score:', as.character(filtered_df()[s, 'MutScore']), 'ADA Score:', as.character(filtered_df()[s, 'ada_score']))
      HTML(paste(patho_vector, collapse = "<br/>"))
    }
  })
  
  output$info10 = renderPlot({ s = input$contents_rows_selected
  data = data.frame(patho = c('Pathogenic', 'Benign'), value= c(filtered_df()[s, 'N_Pathogenic_pred'], filtered_df()[s, 'N_Benign_pred']))
  data$fraction <- data$value / sum(data$value)
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(data$patho, "\n", data$value)
  
  if (sum(data$value)){
    # pie(t(data)[,1])
    
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=patho)) +
      geom_rect() +
      geom_label( x=3.5, aes(y=labelPosition, label=label), size=5) +
      scale_fill_brewer(palette=14) +
      coord_polar(theta="y") +
      xlim(c(2, 4)) +
      theme_void() +
      theme(legend.position = "none")
    # plot_ly(t(data), type = 'pie') %>%
    #   layout(title = 'Pathogenicity Predictors')
  }})
  
  
  
  observeEvent(input$`close-landing-page`, {
    updateTabsetPanel(session, "page", "body-content")
  })
  
  ####################
  ### SAVE SESSION ###
  ####################
  
  
  myBookmarks <- reactiveValues(urlDF = NULL)
  
  observeEvent(input$bookmarkBtn, {
    session$doBookmark()
  })
  
  if (file.exists("/session/app/bookmarks.rds")) {
    myBookmarks$urlDF <- readRDS("/session/app/bookmarks.rds")
  } else {
    myBookmarks$urlDF <- NULL
  }
  
  session$onSessionEnded(function() {
    tmpUrlDF <- isolate({myBookmarks$urlDF})
    
    if (!is.null(tmpUrlDF)) {
      saveRDS(tmpUrlDF, "/session/app/bookmarks.rds")}
  })
  
  setBookmarkExclude(c("bookmarkBtn", "description", "urlTable_cell_clicked", "urlTable_rows_all", "urlTable_rows_current", "urlTable_rows_selected", "urlTable_search", "urlTable_state", "urlTable_row_last_clicked", "contents_state", "contents_rows_selected", "sidebarmenu"))
  
  onBookmarked(fun=function(url){
    if(!url %in% myBookmarks$urlDF$URL){
      if(is.null(myBookmarks$urlDF)){
        myBookmarks$urlDF <- unique(data.table(Description = input$description, URL = paste0("<a href='", url, "'>", url,"</a>"), Timestamp = Sys.time(),  User = Sys.info()[["user"]]), by="URL")
      } else {
        myBookmarks$urlDF <- unique(rbindlist(list(myBookmarks$urlDF, data.table(Description = input$description, URL = paste0("<a href='", url, "'>", url,"</a>"), Timestamp = Sys.time(),  User = Sys.info()[["user"]]))), by="URL")
      }
    }
  })
  
  output$urlTable <- DT::renderDataTable({
    req(myBookmarks$urlDF)
    myBookmarks$urlDF[User %in% Sys.info()[["user"]]]
  }, escape=FALSE)
  
  
  
  #################
  ### DOWNLOADS ###
  #################
  
  
  ###### FILEDATA DOWNLOAD MODULES #####
  
  callModule(downloadObj, id = "download1", data = filtered_df )
  
  ##### DOWNLOAD ANNOTATED FILES #####
  
  output$downloadvcf <- renderUI({
    req(rv$downloadvcf)
    tagList(
      downloadButton("download3", "Download annotate VCF")
    )
  })
  
  output$download3 <- downloadHandler(
    filename = function(){rv$downloadvcf},
    content = function(file){file.copy(paste0("./", rv$downloadvcf), file)}
  )
  
  
}

##### RUN APP #######

enableBookmarking(store = "server")
# runApp(list(ui = ui, server = server), launch.browser= TRUE)

shinyApp(ui, server)



