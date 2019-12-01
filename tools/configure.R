
vsearch_dir <- file.path('inst', 'vsearch')
swarm_dir <- file.path('inst', 'swarm')

# download binaries
if ((!dir.exists(vsearch_dir)) || (!dir.exists(swarm_dir))) {

  unlink(vsearch_dir, recursive = TRUE)
  unlink(swarm_dir, recursive = TRUE)
  
  tmpdir <- 'tmp'
  dir.create(tmpdir, showWarnings = FALSE)
  dir.create('inst', showWarnings = FALSE)

  platform <- Sys.info()[["sysname"]]

  if (grepl('windows', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    # windows binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-win-x86_64.zip'
    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-win-x86_64.zip'
    unwrap <- unzip

  } else if (grepl('linux', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    # linux binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-linux-x86_64.tar.gz'
    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-linux-x86_64.tar.gz'
    unwrap <- untar

  } else {
    # macos binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-macos-x86_64.tar.gz'
    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-macos-x86_64.tar.gz'
    unwrap <- untar

  }
  
  # download vsearch and move binary to inst/vsearch/bin/vsearch(.exe)
  vsearch_bundle <- file.path(tmpdir, basename(vsearch_bin_url))
  download.file(url = vsearch_bin_url, destfile = vsearch_bundle, mode = 'wb')
  unwrap(vsearch_bundle, exdir = tmpdir)
  invisible(file.rename(list.files(tmpdir, pattern = 'vsearch', full.names = TRUE)[1], vsearch_dir))
  if (grepl('windows', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    dir.create(file.path(vsearch_dir, 'bin'))
    invisible(file.rename(file.path(vsearch_dir, 'vsearch.exe'),
                          file.path(vsearch_dir, 'bin', 'vsearch.exe')))
  }
  
  # download swarm and move binary to inst/swarm/bin/swarm(.exe)
  swarm_bundle  <- file.path(tmpdir, basename(swarm_bin_url))
  download.file(url = swarm_bin_url, destfile = swarm_bundle, mode = 'wb')
  unwrap(swarm_bundle, exdir = tmpdir)
  invisible(file.rename(list.files(tmpdir, pattern = 'swarm', full.names = TRUE)[1], swarm_dir))

  # remove tmpdir
  unlink(tmpdir, recursive = TRUE)
}

