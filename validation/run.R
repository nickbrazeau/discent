#!/usr/bin/env Rscript


library(targets)
#............................................................
# check target framework
# NB, can't have this uncommented in real run or script will not work
#...........................................................
#targets::tar_manifest(fields = commands, script = "validation/_targets.R")
#targets::tar_visnetwork(script = "validation/_targets.R")

#............................................................
# actual run
#...........................................................
here::here()
targets::tar_make(script = "/Users/nbrazeau/Documents/Github/discent/validation/_targets.R",
                  store = "/Users/nbrazeau/Documents/Github/discent/validation")
# targets::tar_make_clustermq(workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint
