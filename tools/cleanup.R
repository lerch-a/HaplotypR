

vsearch_dir <- file.path('inst', 'vsearch')
if (dir.exists(vsearch_dir))
  unlink(vsearch_dir, recursive = TRUE)

swarm_dir <- file.path('inst', 'swarm')
if (dir.exists(swarm_dir))
  unlink(swarm_dir, recursive = TRUE)

if (dir.exists('tmp')) 
  unlink('tmp', recursive = TRUE)
  
