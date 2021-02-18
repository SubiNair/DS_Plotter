library(shiny)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(DT)


#RESET TO FILE LOCATIONS AND FILE NAMES
setwd("xxx")
metadata <- data.frame(read.delim(".xxx"))
abundance_data <- data.frame(read.delim("xxx"))
exp_data <- data.frame(read.delim("xxx"))


#extract shared IDs in the metadata
ab_ids <- intersect(unique(metadata$LabID), unique(abundance_data$LabID))
ex_ids <- intersect(unique(metadata$LabID), unique(exp_data$LabID))

#adjusted concentration
ex_data <- filter(metadata, LabID %in% ex_ids)
ex_data_c <- merge(exp_data, ex_data, by="LabID")

#adjusted abundance
ab_data <- filter(metadata, LabID %in% ab_ids)
ab_data_c <- merge(abundance_data, ab_data, by="LabID")





## App

ui <- fluidPage(
  
  
  sidebarPanel(
    h3("Cytokine and Metabolite Regression"),
    
    h5("Data is log2 transformed for plotting purposes."),
    
    selectInput(
      'dataset',
      "Dataset",
      choices = c("Adjusted Relative Abundance" = 'adat', "Adjusted Concentration" = 'edat')
    ), 

    selectInput(
      'karyo',
      "Regression Variable",
      choices = c("Sex" = "Sex", "Age" = "Age", "Karyotype" = "Karyotype")
    ),
    
    textOutput("text")
    
    
    
  ), ### end of sidebar panel
  
  mainPanel(
    
    fluidRow(column(12, div(dataTableOutput("datatb")))),
    plotOutput('plot')
  )
  
)

server <- function(input, output, session) {

  
  cytokines <- reactive({
    req(input$dataset)
    if(input$dataset == 'adat') {
      cks <- unique(ab_data_c$Analyte)
    }
    else{cks <- unique(ex_data_c$Analyte)}
    
    return(cks)
    
  })
  
  
  #Datatable built here
  #log 2 transformation
  
  
  data_table <- reactive({
    req(input$dataset)
    req(input$karyo)
    
    if(input$dataset == 'adat') {
      selected_data <- ab_data_c
      selected_data$adjusted_relative_abundance <- log2(selected_data$adjusted_relative_abundance)
      dframe <- data.frame(matrix(NA, nrow=length(cytokines()), ncol=3))
      colnames(dframe) <- c("Analyte", "Coefficient", "Intercept")
      
      dframe$Analyte <- cytokines()
      
      
      for(val in 1:length(cytokines())) {
        cy_val <- cytokines()[val]
        poi <- dplyr::filter(selected_data, Analyte %in% cy_val)
        model_data <- lm(paste("adjusted_relative_abundance ~", input$karyo), data = poi)
        dframe[val,"Coefficient"] <- model_data$coefficients[2]
        dframe[val,"Intercept"] <- model_data$coefficients[1]
      }
    }
    else {
      selected_data <- ex_data_c
      selected_data$Adjusted_Concentration <- log2(selected_data$Adjusted_Concentration)
      dframe <- data.frame(matrix(NA, nrow=length(cytokines()), ncol=3))
      colnames(dframe) <- c("Analyte", "Coefficient", "Intercept")
      
      dframe$Analyte <- cytokines()
      
      
      for(val in 1:length(cytokines())) {
        cy_val <- cytokines()[val]
        poi <- dplyr::filter(selected_data, Analyte %in% cy_val)
        model_data <- lm(paste("Adjusted_Concentration ~", input$karyo), data = poi)
        dframe[val,"Coefficient"] <- model_data$coefficients[2]
        dframe[val,"Intercept"] <- model_data$coefficients[1]
      }
    }

    return(dframe)
                            
  })
  
  #determine if a boxplot is needed
  x_var <- reactive({
    req(input$karyo)
    req(input$dataset)
    
    colname_s <- input$karyo
    
    if(input$dataset == 'adat') {
      selected_dataset <- ab_data_c
    }
    
    else {
      selected_dataset <- ex_data_c
    }
    
    if(length(unique(selected_dataset[,colname_s])) > 2) {
      box_plot_val_ <- TRUE
    }else{box_plot_val_ <- FALSE}
    
    return(box_plot_val_)
  })

  # plotter reactive object for both the boxplots and standard regression models
  plotter <- reactive({
    req(input$karyo, input$dataset, input$datatb_rows_selected)

  
  
    #extract the selected column
    s = input$datatb_rows_selected
    data_ <- data_table()[s,1]
    
    output$text <- renderText({ paste("Feature selected:", data_, sep = " ") })
    

    if(input$dataset == 'adat') {
      poi <- filter(ab_data_c, Analyte %in% data_)
    }
    else {
      poi <- filter(ex_data_c, Analyte %in% data_)
    }
    
    
    
    box_plot_val <- x_var()
    
    if(isFALSE(box_plot_val)) {
      
      if(input$dataset == 'adat') {
        
        boxplot(adjusted_relative_abundance ~ get(input$karyo) , data=poi, xlab = input$karyo, ylab = "Adjusted Relative Abundance", main = data_)
      }
      else{
         boxplot(Adjusted_Concentration ~ get(input$karyo) , data=poi, xlab = input$karyo, ylab = "Adjusted Concentration", main = data_)
        }
    }
    
    #in the case of age, we use a standard regression
    #can be adjusted for other variables as well, currently supports Age
    
    else {
      if(input$dataset == 'adat') {
        
        g <- ggplot(poi, aes(x = poi$Age, y = poi$adjusted_relative_abundance)) +
          geom_point() +
          stat_smooth(method = lm)
        g + labs(x = input$karyo, y = 'Adjusted Relative Abundance', title = data_)
      }
      else{
        
        g <- ggplot(poi, aes(x = poi$Age, y = poi$Adjusted_Concentration)) +
          geom_point() +
          stat_smooth(method = lm)
        g + labs(x = input$karyo, y = 'Adjusted Concentration', title = data_)
      }
      
    }
  })
  
  output$plot <- renderPlot({
    
    plotter()
  }
  )


  
  output$datatb <- DT::renderDT(
    datatable(data_table(), selection = 'single', rownames = FALSE, options = list(search = list(regex = TRUE), pageLength = 10, autoWidth = TRUE, columnDefs = list(list(width = '200px', targets = "_all")))), server = TRUE)
}

 


shinyApp(ui, server)