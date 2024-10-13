library(rsconnect)

writeManifest(appDir = "app", 
              appPrimaryDoc = "ui.R")

# deployApp(appDir = "app", 
#           #appPrimaryDoc = "ui.R", 
#           appName = "Hurricane_EDA", 
#           appTitle = "Hurricane EDA", 
#           server = "shinyapps.io",
#           account = "tylerpollard410")

deployApp(appDir = "app", 
          appName = "Hurricane_EDA", 
          appTitle = "Hurricane EDA", 
          server = "shinyapps.io",
          account = "tylerpollard410", 
          forceUpdate = TRUE)
