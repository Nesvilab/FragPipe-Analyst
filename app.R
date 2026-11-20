# Single-file entry point for hosts that expect app.R (e.g. Posit Connect Cloud).
# Shiny ignores global.R, ui.R and server.R when app.R exists, so load them here.
source("global.R")
source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)
