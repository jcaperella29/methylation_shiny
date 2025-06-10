shiny::runApp('/app/app.R', host=Sys.getenv('SHINY_HOST'), port=as.numeric(Sys.getenv('SHINY_PORT')), launch.browser=FALSE)
