
homepage <- function(){
  fluidPage(
    titlePanel("Home"),
    hr())
}

configButton <- function(show=F){
  if(show)
    tagList(
        shinySaveButton('saveConfig', 'Save Configuration', 'Save configuration file', filetype=list(rds=c('RDS'))),
        shinyFilesButton('loadConfig', 'Load Configuration', 'Please select the configuration file', FALSE),
        hr()
        )
  else tagList()
}

setinput <- function(){
  fluidPage(
    titlePanel("Input"),
    hr(),
    configButton(getAppOptionDevel()),
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
             h5("Options:"),
             fluidRow(column(2, checkboxInput("withSampleName", label=NULL, value=FALSE)),
                      column(10, h5(helpText("Use sample names for output filenames.")))
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
    h5("Primer sequence selected:"),
    fluidRow(column(12, tableOutput('primerTableDPM'))),
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

concatreads <- function(){
  fluidPage(
    titlePanel("Concatenate Paired Reads"),
    hr(),
    h5("Options:"),
    fluidRow(column(1, checkboxInput("withTrim", label=NULL, value=FALSE)),
             column(11, h5(helpText("Trim reads.")))
    ),
    fluidRow(column(4, h5(helpText("Trim forward read to:"))),
             column(2, numericInput("trim_F", label=NULL, value=200, min=25, step=1)),
             column(4, h5(helpText("Trim reverse read to:"))),
             column(2, numericInput("trim_R", label=NULL, value=200, min=25, step=1))
    ),
    fluidRow(column(2, actionButton("startConcatReads", "Start ...")),
             column(10)
    ),
    hr(),
    h4("Summary:"),
    fluidRow(column(12, tableOutput('table_concat')))
  )
}


callgenotype <- function(){
  fluidPage(
    titlePanel("Call Genotyps"),
    hr(),
    fluidRow(column(2, h5(helpText("Select Marker:"))),
             column(10, uiOutput("uiMarkerSelectionG"))),
    h5("Options:"),
    fluidRow(column(2, numericInput("minMMrate", label=NULL, value=0.5, min=0.5, step=0.05)),
             column(10, h5(helpText("Minimum mis-match rate to call SNP.")))
    ),
    fluidRow(column(2, numericInput("minOccGen", label=NULL, value=2, min=2, step=1)),
             column(10, h5(helpText("Minimum genotype occurence above min mismatch rate to call SNP.")))
    ),
    fluidRow(
      column(2, actionButton("startCallGenotyp", "Start ...")),
      column(10)
    ),
    hr(),
    h4("Plot Average Base Mis-Match Rate per sample:"),
    helpText("(with respect to the reference sequence)"),
    plotOutput('plotMisMatch')
  )
}


callhaplotype <- function(){
  fluidPage(
    titlePanel("Call Haplotypes"),
    hr(),
    fluidRow(column(2, h5(helpText("Select Marker:"))),
             column(10, uiOutput("uiMarkerSelectionH"))),
    h5("Options:"),
    fluidRow(column(1, checkboxInput("cChimera", label=NULL, value=TRUE)),
             column(11, h5(helpText("Check for chimera.")))
    ),
    fluidRow(column(1, checkboxInput("cIndels", label=NULL, value=TRUE)),
             column(11, h5(helpText("Check for indel.")))
    ),
    fluidRow(column(1, checkboxInput("fSNP", label=NULL, value=TRUE)),
    				 column(11, h5(helpText("Use SNP.")))
    ),
    fluidRow(
      column(4, h5(helpText("Minimum coverage per Haplotype:"))),
      column(2, numericInput("minCoverage", label=NULL, value=3, min=3, step=1)),
      column(4, h5(helpText("Detection limit:"))),
      column(2, numericInput("maxSensitivity", label=NULL, value=1/1000, min=0))
    ),
    fluidRow(
      column(4, h5(helpText("Minimum haplotype occurence:"))),
      column(2, numericInput("minOccHap", label=NULL, value=2, min=1, step=1)),
      column(4, h5(helpText("Minimum Sample coverage:"))),
      column(2, numericInput("minCovSample", label=NULL, value=50, min=3, step=1))
    ),
    fluidRow(
      column(2, actionButton("startCallHaplotype", "Start ...")),
      column(10)
    ),
    hr(),
    h4("Summary:"),
    tableOutput('tableHaplotyps')
  )
}
