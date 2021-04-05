library(shiny)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(DT)
library(reshape2)
library(limma)

#RESET TO FILE LOCATIONS AND FILE NAMES
setwd("/Users/subi/Documents/DownSyndromeData")

#comorbidity data
como_data <-read.csv("./ds_data/P4C_Comorbidity_020921.tsv", sep = "\t")
co_meta <- data.frame(read.delim("./ds_data/P4C_metadata_021921_Costello.txt"))

#adjusted concentration
metadata <- data.frame(read.delim("./ds_data/P4C_3omics_metadata_121120.txt"))
ac_data <- data.frame(read.delim("./ds_data/P4C_MSD_Exp3456_JRS_v5.txt"))

#extract shared IDs in the metadata
ex_ids <- intersect(unique(metadata$LabID), unique(ac_data$LabID)) ### 309
co_ids <- intersect(unique(como_data$RecordID), unique(co_meta$RecordID)) ### 468

#309 available IDs to match (ex_ids)

#create merged datasets that contain all of the information required 
#ac_merged contains the 309 lab ids with all of their info from the metadata table as well
ac_in_metadata <- filter(metadata, LabID %in% ex_ids) ### 309
ac_merged <- merge(ac_data, ac_in_metadata, by='LabID') 
### above contains 16686 LabIDs, includes duplicates/multiple rows for same IDs

#patients (RecordIDs) can have multiple visits. the second visit will have a different lab ID.
#como_merged is the same merge but with the comorbidities data
co_data <- filter(co_meta, RecordID %in% co_ids) ### 554 RecordIDs, includes duplicates
como_merged <- merge(como_data, co_data, by='RecordID') ### 7218 rows

#identify top 3 comorbidities and the unique patientIDs that we care about 
top_three <- tail(names(sort(table(como_merged$Condition))), 3)
patientIDs <- unique(como_merged$RecordID)  #468 

#Cut down the dataset -- there are 309 unique LabIDs that we care about
records <- subset(como_merged, LabID %in% ac_merged$LabID)  ### 3834 total, 305 unique LabIDs

#not shared: "HTP0485A" "HTP0487A" "HTP0488A" "HTP0535A"


# records <- records[complete.cases(records), ]
# records = records[!duplicated(records$RecordID),]
# records <- subset(records, Condition %in% top_three)


# length(unique(records$RecordID)) indicates 305, LabIDs and RecordIDs are 1 to 1

x <- data.frame(matrix(nrow=length(unique(ac_merged$Analyte)), ncol = length(unique(records$RecordID))))
colnames(x) <- unique(records$RecordID)
rownames(x) <- unique(ac_merged$Analyte)

# double for loop to go through the column and rows
# Can look into faster options
for(val in colnames(x)){
  lab_id <- unique(records[records$RecordID == val, "LabID"])
  #print(length(unique(lab_id))) indicates 1 to 1, no RecordIDs with multiple LabIDs (otherwise edit above)
  for(cytokine in rownames(x)){
    id_vals <- ac_merged[ac_merged$LabID == lab_id, ]
    index_cyto <- which(id_vals$Analyte == cytokine, arr.ind = TRUE)
    ac_val <- id_vals$Adjusted_Concentration[index_cyto]
    x[cytokine, val] <- ac_val
  }
} #spot check is correct

#Make sure not to run this more than once after this point
x <- log2(x)

#Manually assign for now 
records$Sex[records$Sex == "Male"] <- 1
records$Sex[records$Sex == "Female"] <- 0
records$Karyotype[records$Karyotype == "T21"] <- 1
records$Karyotype[records$Karyotype == "Control"] <- 0
records$Sample_source[records$Sample_source == "NDSC2019"] <- 1
records$Sample_source[records$Sample_source == "Local"] <- 0
records$Event_name[records$Event_name == "Visit 2"] <- 1
records$Event_name[records$Event_name == "Visit 1"] <- 0

#nonnumeric <- subset(records %>% select(1,4,7))
#records <- subset(records %>% select(1:3, 5:6, 8:10))

reg_var <- append(colnames(records)[c(7,8,11,12)], top_three)


## App

ui <- fluidPage(
  
  
  sidebarPanel(
    h3("Cytokine Regression"),
    
    h5("Data is log2 transformed for plotting purposes."),
    
    selectInput(
      'reg',
      "Regression Variable",
      choices = reg_var,
      options = list(maxItems = 1)
    ),
    h5("Two covariates may be selected. Sample source will always be treated as a third covariate."),
    selectizeInput(
      'covar',
      "Covariates",
      choices = NULL,
      multiple = T,
      options = list(maxItems = 2)
    )
    
    
  ), ### end of sidebar panel
  
  mainPanel(
    
    fluidRow(column(12, div(dataTableOutput("datatb")))),
    plotOutput('plot')
  )
  
)

server <- function(input, output, session) {
  
  #Updates the covariates drop down to prevent user from selecting same variable twice
  observeEvent(input$reg, {
    covar_select <- reg_var[!(reg_var %in% input$reg)]

    updateSelectizeInput(session, 'covar',
                         choices =covar_select,
                         server = TRUE,
                         selected = NULL)
  }, once = FALSE)

  
  
  #Data construction
  
  cov_matrix <- reactive({
    req(input$reg)
    req(input$covar)
    
    reg_ <- input$reg
    covar <- input$covar
    
    unprepped <- records[,c('RecordID', 'Condition', 'HasCondition', 'BMI')]
    dc_up <- dcast(unprepped, RecordID + BMI ~ Condition, value.var = 'HasCondition')
    
    # RecordID      BMI Anxiety Any autoimmune skin condition Any congenital heart defect
    # 1 INVAB226VEU 30.81652      NA                            NA                           1
    # 2 INVAC526BJC 16.30712       0                             1                           0
    # 3 INVAE553ERT 14.80047       0                            NA                           1
    # 4 INVAE957TK0 18.14209       0                             0                           0
    # 5 INVAG784LWV 32.76571       0                            NA                           0
    # 6 INVAG895GJN 43.87061       1                             0                           1
    
    dc_up <- subset(dc_up, select=reg_var)
    
  })

  data_table <- reactive({
    req(input$reg)
    req(input$covar)
    
    if(length(setdiff(input$covar, top_three)) > 0) {
      
    }
    
    #first batch correction
    batch1 <- as.factor(ssource)
    #covs <- cbind(A=c(1,1,1,0,0,0), B=c(1,0,1,0,1,0), C=c(20.3,30.4,23.5,32.1,28.5,25))
    #colnames(records) <- c("A1","A2","A3", "A4", "B1","B2","B3", "B4")

    records <- ac_data
    rec_updated <- removeBatchEffect(x, batch = records)
    
    # for(value in input$covar) {
    #   new_batch <- records$value
    #   rec_updated < removeBatchEffect(rec_updated, batch = new_batch)
    # }
    # 
    
    return(rec_updated)


  })
  
  #determine if a boxplot is needed
  x_var <- reactive({
    req(input$reg)
    req(input$covar)
    
    colname_s <- input$reg
    
        if(length(unique(records[,colname_s])) > 2) {
      box_plot_val_ <- TRUE
    }else{box_plot_val_ <- FALSE}
    
    return(box_plot_val_)
  })


  
  output$datatb <- DT::renderDT(
    datatable(data_table(), selection = 'single', rownames = FALSE, options = list(search = list(regex = TRUE), pageLength = 10, autoWidth = TRUE, columnDefs = list(list(width = '200px', targets = "_all")))), server = TRUE)
}

 


shinyApp(ui, server)