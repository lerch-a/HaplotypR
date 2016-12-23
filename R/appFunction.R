
runShinyApp <- function(){
  require(shiny)
  runApp(system.file(package = "haplotypR", "shinyApp", "app.R"))
}
