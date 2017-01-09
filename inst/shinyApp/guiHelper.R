
homepage <- function(){
  fluidPage(
    titlePanel("Home"),
    hr())
}

setinput <- function(){
  fluidPage(
    titlePanel("Input"),
    hr(),
    shinySaveButton('saveConfig', 'Save Configuration', 'Save configuration file', filetype=list(rds=c('RDS'))),
    shinyFilesButton('loadConfig', 'Load Configuration', 'Please select the configuration file', FALSE),
    hr(),
    fluidRow(
      column(12,
             h5("Choose sequence run files (fastq format):"),
             fluidRow(
               column(2, shinyFilesButton('fileR1', 'Browse...', 'Please select the forward read file', FALSE)),
               column(4, textOutput('filepathsR1')),
               column(2, shinyFilesButton('fileR2', 'Browse...', 'Please select the reverse read file', FALSE)),
               column(4, textOutput('filepathsR2'))
             ),
             h5("Choose barcode sequence files (fasta format):"),
             fluidRow(
               column(2, shinyFilesButton('fileB1', 'Browse...', 'Please select the forward barcode file', FALSE)),
               column(4, textOutput('filepathsB1')),
               column(2, shinyFilesButton('fileB2', 'Browse...', 'Please select the reverse barcode file', FALSE)),
               column(4, textOutput('filepathsB2'))
             ),
             # h5("Add linker sequence:"),
             # fluidRow(
             #   column(6, textInput(inputId='linkerF', label=NULL, value="", placeholder="Input forward linker sequence")),
             #   column(6, textInput(inputId='linkerR', label=NULL, value="", placeholder="Input reverse linker sequence"))
             # ),
             h5("Choose primer sequence files (tab delimited format):"),
             fluidRow(
               column(2, shinyFilesButton('fileP1', 'Browse...', 'Please select primer file', FALSE)),
               column(10, textOutput('filepathsP1'))
             ),
             h5("Choose sample files (tab delimited format):"),
             fluidRow(
               column(2, shinyFilesButton('fileSample', 'Browse...', 'Please select sample file', FALSE)),
               column(4, textOutput('filepathsSample')),
               column(6, shinySaveButton('createSampleFile', 'Create Template Sample File', 'Save Sample File', filetype=list(txt=c('txt'))))
             ),
             h5("Choose output directory:"),
             fluidRow(
               column(2, shinyDirButton('dirOut', 'Browse...', 'Please select output directory', FALSE)),
               column(10, textOutput('dirpathOut'))
             )
             # ,
             # fluidRow(
             #   column(2, textInput('iTest', "Test", value="TestA")),
             #   column(10, textOutput('oTest'))
             # ),
             ,hr(),
             fluidRow(uiOutput("inputDyn"))
               
      ))
  )
}

viewSampleList <- function(){
  fluidPage(
    titlePanel("Sample List"),
    hr(),
    fluidRow(column(12, tableOutput('sampleTable')))
    )
}

viewMarkerList <- function(){
  fluidPage(
    titlePanel("Marker List"),
    hr(),
    fluidRow(column(12, tableOutput('markerTable')))
  )
}

demultiplexSample <- function(){
  fluidPage(
    titlePanel("De-multiplex by Sample"),
    hr(),
    fluidRow(
      column(12,
             h5("Sequence run files:"),
             fluidRow(
               column(6, textOutput('fnR1_deplex')),
               column(6, textOutput('fnR2_deplex'))
             ),
             h5("Barcode sequence files:"),
             fluidRow(
               column(6, textOutput('fnB1_deplex')),
               column(6, textOutput('fnB2_deplex'))
             ),
             h5("Sample file:"),
             fluidRow(
               column(12, textOutput('fnS_deplex'))
             ),
             # h5("Choose output directory:"),
             # fluidRow(
             #   column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
             #   column(10, textOutput("dirpathOutDePlex"))
             # ),
             fluidRow(column(2, checkboxInput("withSampleName", label=NULL, value=FALSE)),
                      column(10, h5(helpText("Use sample names for output filenames.")))
             ),
             fluidRow(
               column(2, p(), actionButton("startDePlex", "Start ...")),
               column(10)
             )
      )),
    hr(),
    h4("Summary:"),
    fluidRow(column(12, tableOutput('table_dePlexS')))
  )
}

demultiplexMarker <- function(){
  fluidPage(
    titlePanel("De-Multiplex by Marker"),
    hr(),
    # h5("Linker sequence:"),
    # fluidRow(
    #   column(6, textOutput('linkerF_rmprim')),
    #   column(6, textOutput('linkerR_rmprim'))
    # ),
    h5("Primer sequence selected:"),
    fluidRow(column(12, tableOutput('primerTableDPM'))),
    # h5("Output directory:"),
    # fluidRow(
    #   #column(2, shinyDirButton('dirOutTruncPrim', 'Browse...', 'Please select the output directory', FALSE)),
    #   column(6, textOutput("dirOutTruncPrim")),
    #   column(6)
    # ),
    h5("Options:"),
    fluidRow(
      column(2, numericInput("numMisMatch", label=NULL, value=2, min=0, step=1)),
      column(10, h5(helpText("Number of permited mis-match in primer sequence")))
    ),
    fluidRow(column(2, checkboxInput("withIndel", label=NULL, value=FALSE)),
             column(10, h5(helpText("Permited indels as mis-match in primer sequence")))
    ),
    fluidRow(
      column(2, actionButton("startTruncPrim", "Start ...")),
      column(10)
    ),
    hr(),
    h4("Summary:"),
    fluidRow(column(12, tableOutput('table_dePlexM')))
  )
}
