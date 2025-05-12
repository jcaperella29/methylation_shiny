options(shinyServerMinVersion = "1.0.0")
shiny::runApp("/app", host = Sys.getenv("SHINY_HOST"), port = as.numeric(Sys.getenv("SHINY_PORT")))
