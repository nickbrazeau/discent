#!/usr/bin/env Rscript


library(targets)
#............................................................
# check target framework
#...........................................................
targets::tar_manifest(fields = commands, script = "validation/_targets.R")
targets::tar_visnetwork(script = "validation/_targets.R")

#............................................................
# actual run
#...........................................................
targets::tar_make(script = "validation/_targets.R",
                  store = "validation")
# targets::tar_make_clustermq(workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint
