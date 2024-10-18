#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(bs4Dash)
library(shinyWidgets)
library(waiter)
library(data.table)
library(glue)
library(scales)
library(sf)
library(spData)
library(bestNormalize)
library(smplot2)
library(reactable)
library(lubridate)
library(tidyverse)

# Define UI for application that draws a histogram
shinyUI(
    bs4DashPage(
        title = "Hurricane Explore", skin = "dark", 
        #preloader = list(html = tagList(spin_1(), "Loading ..."), color = "#3c8dbc"),
        
        # Header ====
        header = dashboardHeader(
            title = dashboardBrand(
                title = div(strong("Hurricane EDA"), style = "align:center"),
                color = "primary"
                #image = "https://github.com/TylerPollard410/TylerPollard410.github.io/blob/master/assets/images/profilepic2.jpeg"
            ),
            rightUi = tagList(
                dropdownMenu(
                    #badgeStatus = "info",
                    type = "notifications",
                    icon = icon("circle-info"),
                    headerText = "",
                    notificationItem(
                        inputId = "triggerAuthor",
                        text = "Author: Tyler Pollard",
                        icon =  icon("user"),
                        status = "danger"
                    )
                )
            )
        ),
        
        # Sidebar ====
        sidebar = bs4DashSidebar(
            id = "mysidebar", 
            sidebarMenu(
                id = "mysidebarmenu",
                menuItem(text = "Data", tabName = "data_tab", icon = icon("table")),
                menuItem(text = "Plots", tabName = "plot_tab", icon = icon("chart-simple"))
            )
        ),
        
        # Body ====
        body = dashboardBody(
            tabItems(
                ## Data Tab ----
                tabItem(tabName = "data_tab",
                        fluidPage(
                            h1("Data"),
                            hr(),
                            reactableOutput(outputId = "dataTableOut")
                            
                        ) # end fluidPage
                ), # end plot tab
                
                ## Plot ----
                tabItem(tabName = "plot_tab",
                        fluidPage(
                            fluidRow(
                                column(
                                    width = 3,
                                    box(title = "Data Inputs", 
                                        width = 12,
                                        prettyRadioButtons(
                                            inputId = "plot_type", 
                                            label = "Please select plot type", 
                                            choices = c(
                                                "Scatter" = "scatter_plot",
                                                "Histogram" = "histogram_plot",
                                                #"Boxplot" = "box_plot",
                                                "Map" = "map_plot"
                                            ), 
                                            selected = "scatter_plot", 
                                            status = "info", 
                                            inline = FALSE
                                        ),
                                        hr(),
                                        uiOutput(outputId = "plotDataInputs")
                                    )
                                ),
                                column(
                                    width = 9,
                                    box(title = "Plot Inputs",
                                        width = 12,
                                        dropdownMenu = dropdownButton(
                                            inputId = "plotOptionsButton",
                                            label = "Plot Options",
                                            inline = TRUE,
                                            right = TRUE,
                                            circle = TRUE,
                                            status = "info",
                                            icon = icon("wrench"),
                                            sliderInput(
                                                inputId = "plot_height", 
                                                label = "Plot height",
                                                min = 400, 
                                                max = 1000,
                                                value = 500)
                                        ),
                                        uiOutput(outputId = "OutPlotUI")
                                    ) # end box
                                ) # end column
                            ) # end fluidRow
                        ) # end fluidPage
                ) # end plot tab
            )
        ) # end body
    ) # end dashbaord Page
) # end ShinyUI





