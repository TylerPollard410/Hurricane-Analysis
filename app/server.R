#### Shiny Server Logic
## Tyler Pollard
## Date: 23 Aug 2024



# Read in data ----
## Clean data ----
Stormdata_raw <- fread("_data/E2_data.csv")
Stormdata <- Stormdata_raw |>
    mutate(
        StormID = factor(StormID),
        basin = factor(basin),
        Date = as_datetime(Date, tz = "UTC")
        #Date = as_datetime(format(Date, "%Y-%m-%d %H:%M:%S"),  tz = "UTC")
    )
#as.POSIXct(Stormdata$Date[1])
# str(Stormdata)

## Create training and test data sets ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

## Transform data ----
# Remove not varying 
StormdataTrain2 <- StormdataTrain |>
    select(-lead_time)

# Create Date vars 
dataYears <- year(StormdataTrain2$Date)
dataMonths <- month(StormdataTrain2$Date, label = TRUE)
dataDays <- day(StormdataTrain2$Date)

StormdataTrain3 <- StormdataTrain2 |>
    mutate(
        Year = factor(dataYears, ordered = TRUE),
        Month = dataMonths
    ) |>
    group_by(StormID) |>
    mutate(
        StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours"))
    ) |>
    select(
        StormID,
        Date,
        Year,
        Month,
        StormElapsedTime,
        everything()
    ) |>
    ungroup()

# Adjust Time for table
StormTableData <- Stormdata |>
    mutate(
        Date = Date + hours(4)
    )

# Plot functions ----
## Column type ----
fact_columns <- colnames(StormdataTrain3)[sapply(StormdataTrain3, function(x){is.factor(x)})]
num_columns <- colnames(StormdataTrain3)[sapply(StormdataTrain3, function(x){is.numeric(x)})]



PlotScatter <- function(x, y = "VMAX", color = NULL, fit_line = FALSE, facet = NULL){
    # Color
    if(is.null(color)){
        g_base <- ggplot(data = StormdataTrain3, aes(x = !!sym(x), y = !!sym(y)))
    }else{
        if(is.numeric(StormdataTrain3 |> pull(color))){
            g_base <- ggplot(data = StormdataTrain3, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color), group = !!sym(color))) +
                scale_color_continuous(low = "red", high = "green")
            
        }else{
            g_base <- ggplot(data = StormdataTrain3, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color), group = !!sym(color))) 
        }
    }
    
    # Plot points
    g_point <- g_base + geom_point() + theme_bw()
    
    # Fit Line
    if(fit_line){
        g_fit_line <- g_point +
            sm_statCorr(
                corr_method = "spearman", 
                label_x = ifelse(x == "Date", 
                                 max(StormdataTrain3 |> pull(x)) - months(3),
                                 0.9*max(StormdataTrain3 |> pull(x))),
                legends = TRUE
            )
    }else{
        g_fit_line <- g_point
    }
    
    # Facet Plot
    if(is.null(facet)){
        g_facet <- g_fit_line
    }else{
        g_facet <- g_fit_line +
            facet_wrap(vars(!!sym(facet)))
    }
    
    # Date scale
    if(x == "Date"){
        g_date <- g_facet +
            scale_x_datetime(date_breaks = "month", 
                             date_minor_breaks = "day", 
                             date_labels = "%b-%Y") +
            theme(
                axis.text.x = element_text(angle = 90)
            )
    }else{
        g_date <- g_facet
    }
    
    # Final Plot
    finalPlot <- g_date
    
    finalPlot
}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    # Data Tab ====
    output$dataTableOut <- renderReactable({
        updateData <- StormdataTrain3 |>
            mutate(
                Date = Date + hours(4)
            )
        reactable(data = updateData, 
                  wrap = FALSE, 
                  resizable = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  columns = list(
                      #"StormID",
                      Date = colDef(name = "Date", 
                                    format = colFormat(datetime = TRUE, hour12 = FALSE),
                                    minWidth = 200),
                      StormElapsedTime = colDef(name = "StormElapsedTime", 
                                                #format = colFormat(time = TRUE),
                                                minWidth = 200),
                      #"basin",
                      "LAT" = colDef(name = "LAT", format = colFormat(digits = 4)),
                      "LON" = colDef(name = "LON", format = colFormat(digits = 4)),
                      "MINSLP" = colDef(name = "MINSLP", format = colFormat(digits = 2)),
                      "SHR_MAG" = colDef(name = "SHR_MAG", format = colFormat(digits = 2)),
                      "STM_SPD" = colDef(name = "STM_SPD", format = colFormat(digits = 2)),
                      "SST" = colDef(name = "SST", format = colFormat(digits = 2)),
                      "RHLO" = colDef(name = "RHLO", format = colFormat(digits = 2)),
                      "CAPE1" = colDef(name = "CAPE1", format = colFormat(digits = 2)),
                      "CAPE3" = colDef(name = "CAPE3", format = colFormat(digits = 2)),
                      "SHTFL2" = colDef(name = "SHTFL2", format = colFormat(digits = 2)),
                      "TCOND7002" = colDef(name = "TCOND7002", format = colFormat(digits = 2), minWidth = 150),
                      "INST2" = colDef(name = "INST2", format = colFormat(digits = 2)),
                      "CP1" = colDef(name = "CP1", format = colFormat(digits = 2)),
                      "TCONDSYM2" = colDef(name = "TCONDSYM2", format = colFormat(digits = 2), minWidth = 150),
                      "COUPLSYM3" = colDef(name = "COUPLSYM3", format = colFormat(digits = 2), minWidth = 150),
                      "HWFI" = colDef(name = "HWFI", format = colFormat(digits = 2)),
                      "VMAX_OP_T0" = colDef(name = "VMAX_OP_T0", format = colFormat(digits = 2), minWidth = 150),
                      "HWRF" = colDef(name = "HWRF", format = colFormat(digits = 2)),
                      "VMAX" = colDef(name = "VMAX", format = colFormat(digits = 2))
                  ) # end columns
        )
    })
    
    # Plot Tab ====
    ## Data Input ----
    output$plotDataInputs <- renderUI({
        ### Scatter Plot ----
        if(input$plot_type == "scatter_plot"){
            tagList(
                pickerInput(
                    inputId = "scatter_x",
                    label = "Select x variable",
                    choices = c("Date", num_columns),
                    selected = NULL
                ),
                hr(),
                pickerInput(
                    inputId = "scatter_y",
                    label = "Select y variable",
                    choices = c(num_columns),
                    selected = "VMAX"
                ),
                hr(),
                materialSwitch(
                    inputId = "scatter_color_switch",
                    label = "Color plot by variable?", 
                    value = FALSE, 
                    status = "info"
                ),
                conditionalPanel(condition = 'input.scatter_color_switch',
                                 pickerInput(
                                     inputId = "scatter_color",
                                     label = "Select color variable",
                                     choices = colnames(StormdataTrain3),
                                     selected = NULL, 
                                     multiple = TRUE,
                                     options = pickerOptions(
                                         maxOptions = 1, 
                                         virtualScroll = 600,
                                         dropupAuto = FALSE
                                     )
                                 )
                ),
                hr(),
                materialSwitch(
                    inputId = "scatter_facet_switch",
                    label = "Facet plot by variable?", 
                    value = FALSE, 
                    status = "info"
                ),
                conditionalPanel(condition = 'input.scatter_facet_switch',
                                 pickerInput(
                                     inputId = "scatter_facet",
                                     label = "Select facet variable",
                                     choices = colnames(StormdataTrain3),
                                     selected = NULL, 
                                     multiple = TRUE,
                                     options = pickerOptions(
                                         maxOptions = 1,
                                         virtualScroll = 600,
                                         dropupAuto = FALSE
                                     )
                                 )
                ),
                hr(),
                materialSwitch(
                    inputId = "scatter_fit_line",
                    label = "Fit line?", 
                    value = FALSE, 
                    status = "info"
                )
            )
        }
    })
    
    
    ## Plot Output -----
    ### Box Side bar ----
    
    ### Plot ----
    output$OutPlot <- renderPlot({
        req(input$scatter_y)
        if(input$plot_type == "scatter_plot"){
            scatterPlot <- PlotScatter(
                x = input$scatter_x, 
                y = input$scatter_y, 
                color = if(input$scatter_color_switch){input$scatter_color},
                fit_line = input$scatter_fit_line, 
                facet = if(input$scatter_facet_switch){input$scatter_facet}
            )
        }
        
        plotOut <- scatterPlot
        return(plotOut)
    })
    
    output$OutPlotUI <- renderUI({
        plotOutput(outputId = "OutPlot",
                   width = "100%",
                   height = input$plot_height)
    })
    
})
