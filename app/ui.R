ui <- function(request){
  
  dashboardPagePlus( skin="green",
                     
                     ##############
                     ### HEADER ###
                     ##############
                     
                     header= dashboardHeaderPlus(  
                       title = "PriorR 2.1",
                       left_menu = tagList (
                         dropdownBlock(
                           id="Help",
                           title="",
                           icon="question-circle",
                           actionButton("pdf", "Help", onclick = "window.open('manual.pdf')", style = "width:100px")))),
                     
                     ###############
                     ### SIDEBAR ###
                     ###############
                     
                     dashboardSidebar(disable = FALSE,
                                      width =250,
                                      
                                      ########## SNVs ##########
                                      
                                      # INPUT #
                                      
                                      div(id='tab1_sidebar',
                                          fileInput("file1", "Upload your SNV File",
                                                    multiple = FALSE,
                                                    accept = c("text/csv",
                                                               "text/comma-separated-values,text/plain",
                                                               ".csv")),
                                          
                                          # SIDEBAR #
                                          
                                          sidebarMenu(id = "sidebarmenu",
                                                      
                                                      # NAME #
                                                      useShinyalert(),  # Set up shinyalert
                                                      verbatimTextOutput("name", placeholder = FALSE),
                                                      actionButton("preview", "  Sing in", icon=icon("user"), width = '200px', style="color: #fff; background-color: #9bcd9b; border-color: #f0f0f0 ; font-size: 17px"),
                                                      
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
                                                               
                                                               pickerInput("Biotype", "Biotype:",
                                                                           choices = NULL,
                                                                           selected = NULL,
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE),
                                                               
                                                               pickerInput("Consequence","Consequence:",
                                                                           choices = NULL,
                                                                           selected = NULL,
                                                                           options = list(`actions-box` = TRUE),
                                                                           multiple = TRUE ),
                                                               
                                                               
                                                               fileInput("file2", "Regions",
                                                                         multiple = FALSE,
                                                                         accept = c("text/csv",
                                                                                    "text/comma-separated-values,text/plain",
                                                                                    ".csv"))),
                                                      
                                                      
                                                      # INHERITANCE #
                                                      
                                                      menuItem(
                                                        "Inheritance",           
                                                        tabName = "Inheritance",
                                                        icon=icon("child"),
                                                        selectInput("Inheritance Pattern", "Inheritance Pattern:",
                                                                    c("Autosomal dominant",
                                                                      "Autosomal recessive",
                                                                      "X-linked dominant",
                                                                      "X-linked recesive"))),
                                                      
                                                      # OPTIONS #
                                                      
                                                      menuItem(
                                                        "Options",       
                                                        tabName = "Options",
                                                        icon=icon("list-alt"),
                                                        prettyCheckbox(inputId= "CANONICAL", label = "CANONICAL", value = TRUE, 
                                                                       outline= TRUE,  bigger = TRUE, status = 'success', width = NULL),
                                                        
                                                        prettyCheckbox(inputId="predicted", label = "EXCLUDE PREDICTED", value = FALSE, 
                                                                       outline= TRUE,  status = 'success', width = NULL)
                                                        
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
                                                                             ".csv")),
                                                        fileInput("file4", "Exclude Virtual Panel",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv")),
                                                        selectizeInput("orphanet", "Gene Panel Orphanet",
                                                                       multiple=TRUE,
                                                                       selected = NULL, 
                                                                       options =list(placeholder='Please select a panel'),
                                                                       before_last_dot(list.files(Paneles_Orphanet, pattern='_panel.txt$')))),
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
                                                      
                                                      # LOH #
                                                      
                                                      menuItem(
                                                        "Family info / LOH",
                                                        tabName = "LOH",
                                                        icon=icon("home"),
                                                        textInput("text2", "Probandus", width = '100%'),
                                                        fileInput("file5", "Pedegree File",
                                                                  multiple = FALSE,
                                                                  accept = c("text/csv",
                                                                             "text/comma-separated-values,text/plain",
                                                                             ".csv")),
                                                        prettyCheckbox(inputId="LOH", label = "LOH", value = FALSE, outline= TRUE,  status = 'success', width = NULL)
                                                      ),
                                                      
                                                      
                                                      
                                                      # COMMERCIAL PIPELINE  #
                                                      
                                                      menuItem(
                                                        "Commercial Pipeline",       
                                                        tabName = "Commercial Pipeline",
                                                        icon=icon("project-diagram"),
                                                        prettyCheckbox(inputId= "not-commercial", label = "COLOURING", value = F, 
                                                                       outline= TRUE,  bigger = TRUE, status = 'success', width = NULL),
                                                        
                                                        prettyCheckbox(inputId="not-commercial-filtering", label = "FILTERING", value = F, 
                                                                       outline= TRUE,  status = 'success', width = NULL)
                                                        
                                                      ),
                                                      
                                                      
                                                      # CLASSIFICATION AND ANNOTATION #
                                                      
                                                      menuItem(
                                                        "Classification", 
                                                        icon=icon("th-list"),
                                                        prettyCheckbox(inputId="Classification", label = "Annotate Classified Variants", value = FALSE, 
                                                                       outline= TRUE, fill = TRUE, status = 'success', width = NULL)),
                                                      
                                                      # FREQUENCIES #
                                                      
                                                      menuItem(
                                                        "Frequecies",
                                                        setSliderColor(c("ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen", "ForestGreen"), c(1, 2, 3, 4, 5)),
                                                        icon=icon("users"),
                                                        sliderInput("freq1", "1000G AF", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq2", "ExAC Adjusted AF", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq3", "GnomAD Exomes AF", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq5", "GnomADG AF POPMAX", min = 0, max = 1, value = 0.01),
                                                        sliderInput("freq4", "Maximum AF", min = 0, max = 1, value = 0.01)),
                                                      
                                                      
                                                      menuItem(
                                                        "Pathogenicity",
                                                        setSliderColor("ForestGreen", 1),
                                                        icon=icon("file-medical-alt"),
                                                        sliderTextInput("patho1", "CADD PHRED", choices = seq(from = 75, to = 0, by = -5), selected = 15)),
                                                      
                                                      
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
                                                                 ".csv")),
                                            textInput("text5", "Sample"),
                                            sidebarMenu(id = "sidebarmenu2",
                                                        
                                                        # COLUMNS #
                                                        
                                                        menuItem("Columns_cnv",
                                                                 tabName = "Columns_cnv",
                                                                 icon=icon("columns"),
                                                                 pickerInput("Columns_cnv",   
                                                                             choices = NULL,
                                                                             selected = NULL,
                                                                             options = list(`actions-box` = TRUE),
                                                                             multiple = TRUE )),
                                                        
                                                        
                                                        # OPTIONS #
                                                        
                                                        
                                                        
                                                        # GENE ANNOTATION #
                                                        
                                                        menuItem(
                                                          "Genes",
                                                          tabName = "genes",
                                                          icon=icon("dna"),
                                                          textInput("text", "Your genes"),
                                                          fileInput("file7", "Your Virtual Panel",
                                                                    multiple = FALSE,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain",
                                                                               ".csv")),
                                                          fileInput("file8", "Exclude Virtual Panel",
                                                                    multiple = FALSE,
                                                                    accept = c("text/csv",
                                                                               "text/comma-separated-values,text/plain",
                                                                               ".csv")),
                                                          selectizeInput("orphanet", "Gene Panel Orphanet",
                                                                         multiple=TRUE,
                                                                         selected = NULL, 
                                                                         options =list(placeholder='Please select a panel'),
                                                                         before_last_dot(list.files(Paneles_Orphanet, pattern='_panel.txt$')))),
                                                        
                                                        # CLASSIFICATION #
                                                        
                                                        menuItem(
                                                          "Classification", 
                                                          icon=icon("th-list"),
                                                          prettyCheckbox(inputId="Classification", label = "Annotate Classified Variants", value = FALSE, 
                                                                         outline= TRUE, fill = TRUE, status = 'success', width = NULL)),
                                                        
                                                        menuItem(
                                                          "Frequencies",
                                                          setSliderColor(c("ForestGreen", "ForestGreen", "ForestGreen"), c(1, 2, 3)),
                                                          icon=icon("users"),
                                                          sliderInput("freq1_cnv", "GnomADSV AF POPMAX", min = 0, max = 1, value = 0.01),
                                                          sliderInput("freq2_cnv", "1000g MAX AF", min = 0, max = 1, value = 0.01),
                                                          sliderInput("freq3_cnv", "IMH AF", min = 0, max = 1, value = 0.01))
                                                        # 
                                                        
                                            )
                                        ))
                     ),
                     
                     ############
                     ### BODY ###
                     ############
                     
                     dashboardBody(
                       useShinyjs(),
                       tabsetPanel(
                         id="navbar",
                         tabPanel(title="SNV",id="tab1",value='tab1_val', DTOutput("contents")),
                         tabPanel(title="CNV",id="tab2",value='tab2_val', DTOutput("contents_cnv"))
                       ),
                       tags$head(tags$style(HTML('
                                                 .main-header .logo {
                                                 background-color: #ad1d28;
                                                 font-family:"Georgia", Times, "Times New Roman", serif;
                                                 font-weight: bold;
                                                 font-size: 52 px;
                                                 }
                                                 .navbar{
                                                 background-color: cyan
                                                 }
                                                 '
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
                       dataTableOutput("hpos_table"),
                       uiOutput("modal"),
                       tags$script("Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
                                   Shiny.onInputChange(variableName, null);
                                   });
                                   ")
                       
                       ))
  
  }
