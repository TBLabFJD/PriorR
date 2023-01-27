server <- function(input, output, session) {
  
  ####################
  ### UPLOAD FILES ###
  ####################
  
  # VARIANTS
  
  
  df <- reactive({
    req(input$file1)
    df <- read.table(input$file1$datapath, fill = TRUE, quote = "", header = TRUE,
                     sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
  })
  
  
  df_cnv <- reactive({
    req(input$file6)
    df_cnv <- read.table(input$file6$datapath, fill = TRUE, header = TRUE,
                         sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
  })
  
  
  # PANEL
  
  panel <- reactive({
    req(input$file3)
    panel <- read.csv(input$file3$datapath,
                      sep = ',', header=FALSE)
  })
  
  excluded_panel <- reactive({
    req(input$file4)
    panel <- read.csv(input$file4$datapath,
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
  
  # DIFF FILE
  
  
  diff <- reactive({
    req(input$file1)
    sample_dna <-strsplit(input$file1$name, "_")[[1]][1]
    diff_file <- paste(sample_route, sample_dna, "_bed_1000padding_comparison.diff.sites_in_files", sep="")
    diff <- read.table(diff_file,
                       fill = TRUE, quote = "", header = TRUE,
                       sep = '\t', na.strings=c("",".","NA"), colClasses = NA)
    
  } )
  
  
  fjd_variants <- reactive({fjd_variants <- diff()[diff()["IN_FILE"] == 1,]})
  
  
  #################
  ### FUNCTIONS ###
  #################
  
  ### USER LOGIN ###
  ######
  observeEvent(input$preview, {
    shinyalert("User Name", type = "input") # Show a modal when the button is pressed
  })
  
  
  output$name <- renderText( paste("User:", input$shinyalert),
                             outputArgs = list(
                             ))
  
  
  observeEvent( input$shinyalert, {
    dir.create(file.path(opt$programdir, input$shinyalert))
    setwd(file.path(opt$programdir, input$shinyalert))
  })
  
  getwd()
  
  observeEvent(df(), {
    req(df())
    updatePickerInput(session, inputId = "Columns", choices = colnames(df()), selected = colnames(df())[c(1:4, 6, 11:18, 21:22, 29:33, 52:66, 69 )])
    updatePickerInput(session, inputId = "Consequence", choices = levels(addNA(df()$Consequence)), selected = levels(addNA(df()$Consequence)) )
    updatePickerInput(session, inputId = "Type", choices = levels(df()$VARIANT_CLASS), selected = levels(df()$VARIANT_CLASS) )
    updatePickerInput(session, inputId = "Biotype", choices = levels(df()$BIOTYPE), selected = levels(df()$BIOTYPE)) 
    #, selected = levels(addNA(df$VEP_BIOTYPE))
  })
  
  
  observeEvent(df_cnv(), {
    req(df_cnv())
    updatePickerInput(session, inputId = "Columns_cnv", choices = colnames(df_cnv()), selected = colnames(df_cnv()))
  })
  
  
  
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  ### FILTERS, CHOICES, FREQUENCIES, PANELS ###
  
  
  # FILTERS AND OPTIONS #
  ######
  
  # FILTERS SNVs #
  
  filtered_df <- reactive( {
    
    DF <- df() %>%
      
      select(input$Columns) %>%
      filter(BIOTYPE %in% input$Biotype )  %>%
      filter(Consequence %in% input$Consequence )  %>%
      filter(Genomic_region %in% input$Region ) %>%
      filter(VARIANT_CLASS  %in% input$Type) %>%
      filter(if (input$CANONICAL == FALSE) CANONICAL == CANONICAL else CANONICAL == 'YES') %>%
      filter(if (input$predicted == FALSE) input$predicted == input$predicted  else !grepl("X", Feature)) %>%
      filter(if (input$text == '') SYMBOL == SYMBOL else SYMBOL %in% list(input$text) ) %>%
      filter(if (is.null(input$file3)) SYMBOL == SYMBOL else SYMBOL %in% panel()[,1] ) %>%
      filter(if (is.null(input$file4)) SYMBOL == SYMBOL else SYMBOL %notin% excluded_panel()[,1] ) %>%
      filter( X1000Gp3_AF <= input$freq1 | is.na(X1000Gp3_AF) ) %>%
      filter( ExAC_Adj_AF <= input$freq2 | is.na(ExAC_Adj_AF) ) %>%
      filter( gnomAD_exomes_AF <= input$freq3 | is.na(gnomAD_exomes_AF))  %>%
      filter( gnomADg_AF_POPMAX <= input$freq5 | is.na(gnomADg_AF_POPMAX))  %>%
      filter( MAX_AF <= input$freq4 | is.na(MAX_AF)) %>%
      filter( CADD_PHRED >= input$patho1 | is.na(CADD_PHRED))
    
    # Orphanet gene sets
    
    if (!is.null(input$orphanet)){
      orphanet_path <- paste(Paneles_Orphanet, input$orphanet, ".txt",sep='/')
      gene_set <- read.table( orphanet_path, header=FALSE, sep='\t' )
      DF <- DF %>% filter(SYMBOL %in% gene_set[,1] )
    }
    else{
      DF$SYMBOL = DF$SYMBOL
    }
    # 
    # Bed file upload
    if (!is.null(input$file2)){
      DF$POS <- as.numeric(as.character(DF$POS))
      DF <- left_join( DF, bed(), c("CHROM" = "V1" )) %>% filter( POS >= V2, POS <= V3 )
    }
    
    # Classified variants annotation
    if (input$Classification == TRUE){
      variants_db <- read.table('//qsfjdfile/2t20/USUARIOS/GENETICA/0/PRIORR/Database/variants_priorr.txt', fill = TRUE, 
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
    if(input$'not-commercial-filtering' == TRUE){
      DF <- semi_join(DF, fjd_variants(), by=c("CHROM", "POS"="POS1"))
    }
    
    
    
    # Variant Classification
    DF$Classification.Button = shinyInput( actionButton, nrow(DF), 'button_', label = "Classify" , onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    
    return(DF)
  })
  
  
  
  # FILTERS CNVs #
  
  
  
  filtered_df_cnv <- reactive( {
    
    DF <- df_cnv() %>%
      
      select(input$Columns_cnv) %>%
      filter(if (input$text5 == '') SAMPLE == SAMPLE else SAMPLE %in% list(input$text5))
    # filter(if (is.null(input$file8)) Gene.name == Gene.name else AnnotSV.type == split & Gene.name %in% panel_cnv()[,1] | AnnotSV.type != split) %>%
    # filter(if (is.null(input$file9)) Gene.name == Gene.name else AnnotSV.type == split & Gene.name %notin% excluded_panel_cnv()[,1] | AnnotSV.type != split) %>%
    # filter(SV.type == strsplit( GD_POPMAX_AF, "_")[[1]][3] & GD_POPMAX_AF <= input$freq6 | SV.type != strsplit( GD_POPMAX_AF, "_")[[1]][3])
    # filter( SV.type == X1000g_event & X1000g_max_AF <= input$freq7 | SV.type != X1000g_event ) 
    # filter(if (SV.type == IMH_SV_type) IMH_AF <= input$freq8 | is.na(IMH_AF) else ...)
    
    
    # FREQUENCIES
    #  GD_SV_type = strsplit( GD_POPMAX_AF, "_")[[1]][3]
    
    
    # Orphanet gene sets
    
    # Gene.name and split
    # 
    # if (!is.null(input$orphanet)){
    #   orphanet_path <- paste('/mnt/genetica/ionut/paneles_de_genes/paneles_orphanet/', input$orphanet, ".txt",sep='')
    #   gene_set <- read.table( orphanet_path, header=FALSE, sep='\t' )
    #   DF <- DF %>% filter(Gene.name %in% gene_set[,1] )
    # }
    # else{
    #   DF$SYMBOL = DF$SYMBOL
    # }
    # 
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
  # MODAL CLASSIFICATION POP_UP SNVs
  
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
  
  # ANNOTATION OF VARIANTS TO SNV DB
  
  observeEvent(input$Enter, {
    selectedRow = as.numeric(strsplit(input$select_button, "_")[[1]][2])
    line = paste( filtered_df()[selectedRow,1],  filtered_df()[selectedRow,2],  filtered_df()[selectedRow,3], filtered_df()[selectedRow,4], input$pathogenicity, input$comments, input$shinyalert, sep='\t')
    print(line)
    write(line, '//qsfjdfile/2t20/USUARIOS/GENETICA/0/PRIORR/Database/variants_priorr.txt', append = TRUE)
    
  })
  
  # LOH #
  
  yellowRows <- reactive({
    
    req(filtered_df())
    
    # Only probandus
    if( input$text2 != ""){
      samples = input$text2
      loh = paste( samples, "LOH", sep= "_")
      print("probandus is:")
      print(loh)
      ind = names(filtered_df()) %in% loh
      which(filtered_df()[, ind] == "True")  - 1L
      
      # With ped file
    }else if (!is.null(input$file5)){
      afected = ped()[ped()$V6 == 2, 2]
      samples1 = lapply(afected, gsub, pattern ="-", replacement=".")
      samples = paste("X", samples1, sep="")
      
      if(length(samples) > 1){
        print("Here I am")
        loh = lapply(samples, paste, "LOH", sep = "_")
        ind = names(filtered_df()) %in% loh
        which(rowSums(filtered_df()[, ind] == "True") == sum(ind))  - 1L
        
      }else{
        loh = paste( samples, "LOH", sep= "_")
        ind = names(filtered_df()) %in% loh
        which(filtered_df()[, ind] == "True")  - 1L}
      
      # Single sample  (only one _LOH column)
    }else {
      sample_loh = names(filtered_df())[endsWith(names(filtered_df()), "_LOH")]
      # print("No ped file needed, the only sample is:")
      # print(sample_loh)
      if(length(sample_loh) == 1 ){
        samples = strsplit(sample_loh, "_")[[1]][1]
        loh = paste( samples, "LOH", sep= "_")
        ind = names(filtered_df()) %in% loh
        which(filtered_df()[, ind] == "True")  - 1L
      }else{
        return(NULL)}
    }
    
  })
  
  print(yellowRows)
  
  
  
  ### TOGGLE MENU ###
  
  ######
  
  values <- reactiveValues(selectedTab = 1)
  
  observeEvent( input$navbar, {
    toggle("tab1_sidebar", condition = input$navbar == "tab1_val")
    toggle("tab2_sidebar", condition = input$navbar == "tab2_val")
  })
  
  
  ### OUTPUT TABLE ###
  
  output$contents <- renderDT({
    req(filtered_df())
    datatable(
      filtered_df(),
      class = "display nowrap compact", 
      filter = "top", 
      callback = JS(callback(yellowRows())),
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = '10px', targets = c(3))),
        pageLength = 10),
      escape = FALSE)},
    server = FALSE)
  
  # Output CNVs 
  
  output$contents_cnv <- renderDT({
    req(filtered_df_cnv())
    datatable2(
      filtered_df_cnv(),
      class = "display nowrap compact", 
      filter = "top", 
      options = list(
        scrollX = TRUE, columnDefs = list(list(autoWidth = TRUE, width = '10px', targets = c(3))),
        pageLength = 10),      
      escape = FALSE) %>% 
      formatStyle( "AnnotSV.ranking", 
                   backgroundColor = styleEqual( c(1,2,3,4,5), c("chartreuse3", "greenyellow", "yellow", "orange", "red" )))
  },
  
  server = FALSE)
  
  ### SAVE SESSION ###
  
  ######
  
  
  myBookmarks <- reactiveValues(urlDF = NULL)
  
  observeEvent(input$bookmarkBtn, {
    session$doBookmark()
  })
  
  if (file.exists("bookmarks.rds")) {
    myBookmarks$urlDF <- readRDS("bookmarks.rds")
  } else {
    myBookmarks$urlDF <- NULL
  }
  
  session$onSessionEnded(function() {
    tmpUrlDF <- isolate({myBookmarks$urlDF})
    
    if (!is.null(tmpUrlDF)) {
      saveRDS(tmpUrlDF, "bookmarks.rds")}
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
  
  ###### FILEDATA DOWNLOAD MODULES #####
  
  callModule(downloadObj, id = "download1", data = filtered_df )
  
  
}