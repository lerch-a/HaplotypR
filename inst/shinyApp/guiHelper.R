
homepage <- function(){
  fluidPage(
    titlePanel("Home"),
    hr())
}

setInput <- function(){
  fluidPage(
    titlePanel("Input"),
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
             h5("Add adapter sequence:"),
             fluidRow(
               column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward adapter sequence")),
               column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse adapter sequence"))
             ),
             h5("Choose primer sequence files (fasta format):"),
             fluidRow(
               column(2, shinyFilesButton('fileP1', 'Browse...', 'Please select the forward primer file', FALSE)),
               column(4, textOutput('filepathsP1')),
               column(2, shinyFilesButton('fileP2', 'Browse...', 'Please select the reverse primer file', FALSE)),
               column(4, textOutput('filepathsP2'))
             ),
             h5("Choose output directory:"),
             fluidRow(
               column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
               column(4, textOutput("dirpathOutDePlex")),
               column(6)
             )

      )),
    hr(),
    fluidRow(column(12, dataTableOutput('sampleTable')))
  )
}


demultiplex <- function(){
  fluidPage(
    titlePanel("De-multiplex Amplicon Sequence Run"),
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
             h5("Add adapter sequence:"),
             fluidRow(
               column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward adapter sequence")),
               column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse adapter sequence"))
             ),
             h5("Choose output directory:"),
             fluidRow(
               column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
               column(4, textOutput("dirpathOutDePlex")),
               column(2, numericInput("numMisMatch", "Mis-match", 0, min=0, step=1)),
               column(4, checkboxInput("withIndel", "with-indel"))
             ),
             fluidRow(
               column(2, actionButton("startDePlex", "Start ...")),
               column(4, verbatimTextOutput("runDePlex")),
               column(2),
               column(4)
             )
      )),
    hr(),
    fluidRow(column(12, dataTableOutput('sampleTable')))
  )
}

rmprimer <- function(){
  fluidPage(
    titlePanel("Remove Primer"),
    hr(),
    fluidRow(
      column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward linker sequence")),
      column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse linker sequence"))
    ),
    h5("Choose primer sequence files (fasta format):"),
    fluidRow(
      column(2, shinyFilesButton('fileP1', 'Browse...', 'Please select the forward primer file', FALSE)),
      column(4, textOutput('filepathsP1')),
      column(2, shinyFilesButton('fileP2', 'Browse...', 'Please select the reverse primer file', FALSE)),
      column(4, textOutput('filepathsP2'))
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
             h5("Add adapter sequence:"),
             fluidRow(
               column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward adapter sequence")),
               column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse adapter sequence"))
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
               column(4, helpText("Select forward marker adapter file")),
               column(2, shinyFilesButton('fileB2', 'Browse...', 'Please select the reverse marker file', FALSE)),
               column(4, helpText("Select reverse marker adapter file"))
             )
      )),
    hr()
  )
}