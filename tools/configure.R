
vsearch_dir <- file.path('inst', 'vsearch')
swarm_dir <- file.path('inst', 'swarm')

# download binaries
if ((!dir.exists(vsearch_dir)) || (!dir.exists(swarm_dir))) {

  unlink(vsearch_dir, recursive = TRUE)
  unlink(swarm_dir, recursive = TRUE)

  tmpdir <- 'tmp'
  dir.create(tmpdir, showWarnings = FALSE)

  platform <- Sys.info()[["sysname"]]

  if (grepl('windows', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    # windows binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-win-x86_64.zip'
    vsearch_bin_zip <- file.path(tmpdir, basename(vsearch_bin_url))
    download.file(url = vsearch_bin_url, destfile = vsearch_bin_zip)
    unzip(vsearch_bin_zip, exdir = tmpdir)

    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-win-x86_64.zip'
    swarm_bin_zip <- file.path(tmpdir, basename(swarm_bin_url))
    download.file(url = swarm_bin_url, destfile = swarm_bin_zip)
    unzip(swarm_bin_zip, exdir = tmpdir)

  } else if (grepl('linux', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    # linux binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-linux-x86_64.tar.gz'
    vsearch_bin_tar <- file.path(tmpdir, basename(vsearch_bin_url))
    download.file(url = vsearch_bin_url, destfile = vsearch_bin_tar)
    untar(vsearch_bin_tar, exdir = tmpdir)

    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-linux-x86_64.tar.gz'
    swarm_bin_tar <- file.path(tmpdir, basename(swarm_bin_url))
    download.file(url = swarm_bin_url, destfile = swarm_bin_tar)
    untar(swarm_bin_tar, exdir = tmpdir)

  } else {
    # macos binaries
    vsearch_bin_url <- 'https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-macos-x86_64.tar.gz'
    vsearch_bin_tar <- file.path(tmpdir, basename(vsearch_bin_url))
    download.file(url = vsearch_bin_url, destfile = vsearch_bin_tar)
    untar(vsearch_bin_tar, exdir = tmpdir)

    swarm_bin_url <- 'https://github.com/torognes/swarm/releases/download/v3.0.0/swarm-3.0.0-macos-x86_64.tar.gz'
    swarm_bin_tar <- file.path(tmpdir, basename(swarm_bin_url))
    download.file(url = swarm_bin_url, destfile = swarm_bin_tar)
    untar(swarm_bin_tar, exdir = tmpdir)

  }

  # move vsearch install to 'inst/vsearch'
  dir.create('inst', showWarnings = FALSE)

  # rename vsearch and include src
  invisible(file.rename(list.files(tmpdir, pattern = 'vsearch', full.names = TRUE)[1], vsearch_dir))
  vsearch_src_url <- 'https://github.com/torognes/vsearch/archive/v2.14.1.tar.gz'
  download.file(url = vsearch_src_url, destfile = file.path(vsearch_dir, basename(vsearch_src_url)))
  if (grepl('windows', Sys.info()[["sysname"]], ignore.case = TRUE)) {
    dir.create(file.path(vsearch_dir, 'bin'))
    invisible(file.rename(file.path(vsearch_dir, 'vsearch.exe'),
                          file.path(vsearch_dir, 'bin', 'vsearch.exe')))
  }

  # rename vsearch and include src
  invisible(file.rename(list.files(tmpdir, pattern = 'swarm', full.names = TRUE)[1], swarm_dir))
  swarm_src_url <- 'https://github.com/torognes/swarm/archive/v3.0.0.tar.gz'
  download.file(url = swarm_src_url, destfile = file.path(swarm_dir, basename(swarm_src_url)))

  # remove tmpdir
  unlink(tmpdir, recursive = TRUE)
}

