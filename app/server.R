#### Shiny Server Logic
## Tyler Pollard
## Date: 23 Aug 2024



# Read in data ----
## Clean data ----
Stormdata_raw <- fread("_data/E2_data.csv")
## Clean data ----
Stormdata1 <- Stormdata_raw |>
    mutate(
        Obs = 1:nrow(Stormdata_raw),
        StormID = factor(StormID),
        basin = factor(basin),
        Date = as_datetime(Date, tz = "UTC"),
        Year = year(Date),
        Month = month(Date, label = TRUE),
        Day = yday(Date),
        LON2 = LON - 360,
        DataType = ifelse(is.na(VMAX), "Test", "Train")
    ) |>
    group_by(StormID) |>
    mutate(
        StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours"))
    ) |>
    ungroup() |>
    select(
        DataType,
        StormID,
        Date,
        Year,
        Month,
        Day,
        basin,
        StormElapsedTime,
        LAT,
        LON,
        LON2,
        everything(),
        -lead_time
    )

# Create Date vars
# dataYears <- year(Stormdata$Date)
# dataMonths <- month(Stormdata$Date, label = TRUE)
# dataDays <- day(Stormdata$Date)
# dataYearDay <- yday(Stormdata$Date)

### Land Indicator ----
pts <- st_as_sf(Stormdata1, # |> select(LON, LAT), 
                coords = c("LON2", "LAT"),
                crs = 4326
)
land_pts <- !is.na(as.numeric(st_intersects(pts, world)))

Stormdata <- Stormdata1 |>
    mutate(
        Land = factor(land_pts, 
                      labels = c("Water", "Land"))
    ) |>
    select(
        DataType,
        StormID,
        Date,
        Year,
        Month,
        Day,
        basin,
        StormElapsedTime,
        LAT,
        LON,
        Land,
        everything(),
        -LON2,
        -Obs
    )

## Create training and test data sets ----
StormdataTrain <- Stormdata |> filter(complete.cases(VMAX))
StormdataTest <- Stormdata |> filter(!complete.cases(VMAX))

## Transform data ----
# Remove not varying 
# StormdataTrain2 <- StormdataTrain |>
#     select(-lead_time)
# 
# # Create Date vars 
# dataYears <- year(StormdataTrain2$Date)
# dataMonths <- month(StormdataTrain2$Date, label = TRUE)
# dataDays <- day(StormdataTrain2$Date)
# dataYearDays <- yday(StormdataTrain2$Date)

StormdataTrain3 <- StormdataTrain |>
    # mutate(
    #     Year = factor(dataYears, ordered = TRUE),
    #     Month = dataMonths,
    #     Day = dataYearDays
    # ) |>
    # group_by(StormID) |>
    # mutate(
    #     StormElapsedTime = as.numeric(difftime(Date, min(Date), units = "hours"))
    # ) |>
    # select(
    #     StormID,
    #     Date,
    #     Year,
    #     Month,
    #     Day,
    #     StormElapsedTime,
    #     everything()
    # ) |>
    # ungroup() |>
    mutate(
        #"log(VMAX)" = log(VMAX),
        "VMAX/HWRF" = VMAX/HWRF,
        #"log(VMAX/HWRF)" = log(VMAX/HWRF)
    )

# Adjust Time for table
StormTableData <- Stormdata |>
    mutate(
        Date = Date + hours(4)
    )

# Plot functions ----
# Column Type
fact_columns <- colnames(StormdataTrain3)[sapply(StormdataTrain3, function(x){is.factor(x)})]
num_columns <- colnames(StormdataTrain3)[sapply(StormdataTrain3, function(x){is.numeric(x)})]

# Scatter Plot
PlotScatter <- function(x, y = "VMAX", transX = "None", scaleX = FALSE, transY = "None", scaleY = FALSE, color = NULL, fit_line = FALSE, facet = NULL){
    # plotData <- StormdataTrain3 |>
    #     rowwise() |>
    #     mutate(
    #         !!sym(x) := ifelse(transX == "None", !!sym(x),
    #                            ifelse(transX == "Log", log(!!sym(x)), scale(!!sym(x)))),
    #         !!sym(y) := ifelse(transY == "None", !!sym(y),
    #                            ifelse(transX == "Log", log(!!sym(y)), scale(!!sym(y))))
    #     ) 
    plotData <- StormdataTrain3
    
    # Transform X
    if(transX == "Log"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){log(x)})
            )
    }else if(transX == "Arcsinh"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){log(x + sqrt(x^2 + 1))})
            )
    }else if(transX == "YeoJohnson"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){predict(yeojohnson(x, standardize = FALSE))})
            )
    }
    
    # Scale X
    if(scaleX){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){scale(x)})
            )
    }
    
    # Transform Y
    if(transY == "Log"){
        plotData <- plotData |>
            mutate(
                across(!!sym(y), function(x){log(x)})
            )
    }else if(transY == "Arcsinh"){
        plotData <- plotData |>
            mutate(
                across(!!sym(y), function(x){log(x + sqrt(x^2 + 1))})
            )
    }else if(transX == "YeoJohnson"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){predict(yeojohnson(x, standardize = FALSE))})
            )
    }
    
    # Scale Y
    if(scaleY){
        plotData <- plotData |>
            mutate(
                across(!!sym(y), function(x){scale(x)})
            )
    }
    
    # Color
    if(is.null(color)){
        g_base <- ggplot(data = plotData, aes(x = !!sym(x), y = !!sym(y)))
    }else{
        if(is.numeric(plotData |> pull(color))){
            g_base <- ggplot(data = plotData, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color), group = !!sym(color))) +
                scale_color_continuous(low = "red", high = "green")
            
        }else{
            g_base <- ggplot(data = plotData, aes(x = !!sym(x), y = !!sym(y), color = !!sym(color), group = !!sym(color))) 
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
                                 max(plotData |> pull(x)) - months(3),
                                 0.9*max(plotData |> pull(x))),
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
            theme_bw() +
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

# Scatter Plot
PlotHistogram <- function(x = "VMAX", transX = "None", scaleX = FALSE, density = FALSE, facet = NULL){
    plotData <- StormdataTrain3
    
    # Transform X
    if(transX == "Log"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){log(x)})
            )
    }else if(transX == "Arcsinh"){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){log(x + sqrt(x^2 + 1))})
            )
    }
    
    # Scale X
    if(scaleX){
        plotData <- plotData |>
            mutate(
                across(!!sym(x), function(x){scale(x)})
            )
    }
    
    # Histogram
    g_hist <- ggplot(data = plotData) +
        geom_histogram(
            aes(x = !!sym(x), after_stat(density)),
            color = "#99c7c7", fill = "#bcdcdc") +
        theme_bw()
    
    # Density
    if(density){
        g_density <- g_hist +
            geom_density(
                aes(x = !!sym(x)),
                color = "#007C7C", 
                linewidth = 1
            )
    }else{
        g_density <- g_hist
    }
    
    # Facet Plot
    if(is.null(facet)){
        g_facet <- g_density
    }else{
        g_facet <- g_density +
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

# Map plot
world_coordinates <- map_data("world") 
mapPlot <- ggplot() + 
    # geom_map() function takes world coordinates  
    # as input to plot world map 
    geom_map( 
        data = world_coordinates, map = world_coordinates, 
        aes(map_id = region) 
    ) 

PlotMap <- function(color = NULL){
    # Color
    if(is.null(color)){
        map_plot1 <- mapPlot +
            geom_point(
                data = StormdataTrain3,
                aes(x = LON-360, y = LAT)) +
            xlim(c(-180,0)) +
            ylim(c(0,60)) 
    }else{
        if(is.numeric(StormdataTrain3 |> pull(color))){
            map_plot1 <- mapPlot +
                geom_point(
                    data = StormdataTrain3,
                    aes(x = LON-360, y = LAT, color = !!sym(color))) +
                xlim(c(-180,0)) +
                ylim(c(0,60)) +
                scale_color_continuous(low = "green", high = "red")
        }else{
            map_plot1 <- mapPlot +
                geom_point(
                    data = StormdataTrain3,
                    aes(x = LON-360, y = LAT, color = !!sym(color))) +
                xlim(c(-180,0)) +
                ylim(c(0,60))
        }
    }
    
    # # Facet 
    # if(is.null(facet)){
    #     map_plot2 <- map_plot1
    # }else{
    #     map_plot2 <- map_plot1 +
    #         facet_wrap(vars(!!sym(facet)))
    # }
    
    map_plot_final <- map_plot1 +
        labs(
            x = "Longitude",
            y = "Latitude"
        )
    theme_bw()
    
    return(map_plot_final)
}

# Define server logic required to draw a histogram
shinyServer(function(input, output, session){
    
    # Data Tab ----
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
    
    # Plot Tab ----
    ## Data Input ----
    
    output$plotDataInputs <- renderUI({
        ### Scatter Plot ----
        if(input$plot_type == "scatter_plot"){
            tagList(
                pickerInput(
                    inputId = "scatter_x",
                    label = "Select x variable",
                    choices = c(num_columns,"Date"),
                    selected = "HWRF"
                ),
                radioGroupButtons(
                    inputId = "scatter_transform_x",
                    label = "Transform X",
                    choices = c("None",
                                "Log", 
                                "Arcsinh",
                                "YeoJohnson"),
                    individual = TRUE,
                    checkIcon = list(
                        yes = tags$i(class = "fa fa-check-square", 
                                     style = "color: steelblue"),
                        no = tags$i(class = "fa fa-circle-o", 
                                    style = "color: steelblue"))
                ),
                materialSwitch(
                    inputId = "scatter_scale_x",
                    label = "Scale X", 
                    value = FALSE, 
                    status = "info"
                ),
                hr(),
                pickerInput(
                    inputId = "scatter_y",
                    label = "Select y variable",
                    choices = c(num_columns),
                    selected = "VMAX"
                ),
                radioGroupButtons(
                    inputId = "scatter_transform_y",
                    label = "Transform Y",
                    choices = c("None",
                                "Log", 
                                "Arcsinh",
                                "YeoJohnson"),
                    individual = TRUE,
                    checkIcon = list(
                        yes = tags$i(class = "fa fa-check-square", 
                                     style = "color: steelblue"),
                        no = tags$i(class = "fa fa-circle-o", 
                                    style = "color: steelblue"))
                ),
                materialSwitch(
                    inputId = "scatter_scale_y",
                    label = "Scale Y", 
                    value = FALSE, 
                    status = "info"
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
        ## Histogram Plot ----
        else if(input$plot_type == "histogram_plot"){
            tagList(
                pickerInput(
                    inputId = "histogram_x",
                    label = "Select x variable",
                    choices = c(num_columns,"Date"),
                    selected = "VMAX"
                ),
                radioGroupButtons(
                    inputId = "histogram_transform_x",
                    label = "Transform X",
                    choices = c("None",
                                "Log", 
                                "Arcsinh"),
                    selected = "None",
                    individual = TRUE,
                    checkIcon = list(
                        yes = tags$i(class = "fa fa-check-square", 
                                     style = "color: steelblue"),
                        no = tags$i(class = "fa fa-circle-o", 
                                    style = "color: steelblue"))
                ),
                materialSwitch(
                    inputId = "histogram_scale_x",
                    label = "Scale X", 
                    value = FALSE, 
                    status = "info"
                ),
                # hr(),
                # materialSwitch(
                #     inputId = "histogram_color_switch",
                #     label = "Color plot by variable?", 
                #     value = FALSE, 
                #     status = "info"
                # ),
                # conditionalPanel(condition = 'input.histogram_color_switch',
                #                  pickerInput(
                #                      inputId = "histogram_color",
                #                      label = "Select color variable",
                #                      choices = colnames(StormdataTrain3),
                #                      selected = NULL, 
                #                      multiple = TRUE,
                #                      options = pickerOptions(
                #                          maxOptions = 1, 
                #                          virtualScroll = 600,
                #                          dropupAuto = FALSE
                #                      )
                #                  )
                # ),
                hr(),
                materialSwitch(
                    inputId = "histogram_facet_switch",
                    label = "Facet plot by variable?", 
                    value = FALSE, 
                    status = "info"
                ),
                conditionalPanel(condition = 'input.histogram_facet_switch',
                                 pickerInput(
                                     inputId = "histogram_facet",
                                     label = "Select facet variable",
                                     choices = c("basin", "Land", "StormID"),
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
                    inputId = "histogram_density",
                    label = "Overlay density?", 
                    value = FALSE, 
                    status = "info"
                )
            )
        }else if(input$plot_type == "map_plot"){
            ## Map Plot ----
            tagList(
                materialSwitch(
                    inputId = "map_color_switch",
                    label = "Color plot by variable?", 
                    value = FALSE, 
                    status = "info"
                ),
                conditionalPanel(condition = 'input.map_color_switch',
                                 pickerInput(
                                     inputId = "map_color",
                                     label = "Select color variable",
                                     choices = colnames(StormdataTrain3 |> select(-LAT, -LON)),
                                     selected = "VMAX", 
                                     multiple = FALSE,
                                     options = pickerOptions(
                                         #maxOptions = 1, 
                                         virtualScroll = 600,
                                         dropupAuto = FALSE
                                     )
                                 )
                )
            )
        }
    })
    
    ## Plot Output ----
    output$OutPlot <- renderPlot({
        req(input$plot_type)
        if(input$plot_type == "scatter_plot"){
            req(
                input$scatter_x,
                input$scatter_y
            )
            scatterPlot <- PlotScatter(
                x = input$scatter_x,
                y = input$scatter_y,
                transX = input$scatter_transform_x,
                scaleX = input$scatter_scale_x,
                transY = input$scatter_transform_y,
                scaleY = input$scatter_scale_y,
                color = if(input$scatter_color_switch){input$scatter_color},
                fit_line = input$scatter_fit_line,
                facet = if(input$scatter_facet_switch){input$scatter_facet}
            )
            plotOut <- scatterPlot
        }else if(input$plot_type == "histogram_plot"){
            req(
                input$histogram_x
            )
            histogramPlot <- PlotHistogram(
                x = input$histogram_x,
                transX = input$histogram_transform_x,
                scaleX = input$histogram_scale_x,
                #color = if(input$histogram_color_switch){input$histogram_color},
                density = input$histogram_density,
                facet = if(input$histogram_facet_switch){input$histogram_facet}
            )
            plotOut <- histogramPlot
        }else if(input$plot_type == "map_plot"){
            req(
                !is.null(input$map_color_switch)
            )
            mapPlot <- PlotMap(
                color = if(input$map_color_switch){input$map_color}
            )
            plotOut <- mapPlot
        }
        
        return(plotOut)
    })
    
    output$OutPlotUI <- renderUI({
        req(input$plot_type)
        plotOutput(outputId = "OutPlot",
                   width = "100%",
                   height = input$plot_height)
    })
    
})
