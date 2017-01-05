
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
             h5("Add linker sequence:"),
             fluidRow(
               column(6, textInput(inputId='linkerF', label=NULL, value="", placeholder="Input forward linker sequence")),
               column(6, textInput(inputId='linkerR', label=NULL, value="", placeholder="Input reverse linker sequence"))
             ),
             h5("Choose primer sequence files (tab delimited format):"),
             fluidRow(
               column(2, shinyFilesButton('fileP1', 'Browse...', 'Please select primer file', FALSE)),
               column(10, textOutput('filepathsP1'))
             ),
             h5("Choose output directory:"),
             fluidRow(
               column(2, shinyDirButton('dirOut', 'Browse...', 'Please select output directory', FALSE)),
               column(10, textOutput('dirpathsOut'))
             )
      ))
  )
}


demultiplex <- function(){
  fluidPage(
    titlePanel("De-multiplex Sample"),
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
             h5("Choose output directory:"),
             fluidRow(
               column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
               column(10, textOutput("dirpathOutDePlex"))
             ),
             fluidRow(
               column(2, p(), actionButton("startDePlex", "Start ...")),
               column(10)
             )
      )),
    hr(),
    fluidRow(column(12, dataTableOutput('sampleTable_deplex')))
  )
}

demultiplexMarker <- function(){
  fluidPage(
    titlePanel("De-multiplex Marker"),
    hr(),
    h5("Linker sequence:"),
    fluidRow(
      column(6, textOutput('linkerF_rmprim')),
      column(6, textOutput('linkerR_rmprim'))
    ),
    h5("Primer sequence files:"),
    fluidRow(
      column(6, textOutput('fnP1_rmprim')),
      column(6, textOutput('fnP2_rmprim'))
    ),
    h5("Choose output directory:"),
    fluidRow(
      column(2, shinyDirButton('dirOutTruncPrim', 'Browse...', 'Please select the output directory', FALSE)),
      column(4, textOutput("dirOutTruncPrim")),
      column(2, numericInput("numMisMatch", "Mis-match", 0, min=0, step=1)),
      column(4, checkboxInput("withIndel", "with-indel"))
    ),
    fluidRow(
      column(2, actionButton("startTruncPrim", "Start ...")),
      column(4, verbatimTextOutput("runTruncPrim")),
      column(2),
      column(4)
    )
  )
}

temp <- function(){
  fluidPage(
    titlePanel(""),
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
             h5("Add linker sequence:"),
             fluidRow(
               column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward linker sequence")),
               column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse linker sequence"))
             ),
             h5("Choose xxx:"),
             fluidRow(
               column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
               column(4, textOutput('dirpathOutDePlex')), 
               column(2, numericInput("numMisMatch", "Mis-match", 0, min=0, step=1)),
               column(4, checkboxInput("withIndel", "with-indel"),
                      actionButton("startDePlex", "Start ..."),
                      verbatimTextOutput("runDePlex"))
             ),
             h5("Choose marker sequence files (fasta format):"),
             fluidRow(
               column(2, shinyFilesButton('fileB1', 'Browse...', 'Please select the forward marker file', FALSE)),
               column(4, helpText("Select forward marker linker file")),
               column(2, shinyFilesButton('fileB2', 'Browse...', 'Please select the reverse marker file', FALSE)),
               column(4, helpText("Select reverse marker linker file"))
             )
      )),
    hr()
  )
}