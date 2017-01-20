
runShinyApp <- function(){
  require(shiny)
  runApp(system.file(package = "HaplotypR", "shinyApp", "app.R"))
}
