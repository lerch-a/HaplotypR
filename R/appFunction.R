
runShinyApp <- function(devel=FALSE){
  require(shiny)
  if(devel)
    options(HaplotypR.devel=TRUE)
  runApp(system.file(package = "HaplotypR", "shinyApp", "app.R"))
}

getAppOptionDevel <- function(){
  opt <- getOption("HaplotypR.devel")
  if(is.null(opt))
    return(FALSE)
  else
    return(opt==TRUE)
}

