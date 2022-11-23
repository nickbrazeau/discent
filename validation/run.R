#!/usr/bin/env Rscript


library(targets)
#............................................................
# check target framework
# NB, can't have this uncommented in real run or script will not work
#...........................................................
#targets::tar_manifest(fields = commands, script = "validation/_targets_sim.R")
#targets::tar_visnetwork(script = "validation/_targets_sim.R")

#............................................................
# actual run
#...........................................................
setwd("/proj/ideel/meshnick/users/NickB/Projects/discent")
targets::tar_make_clustermq(workers = 24,
                            script = "validation/_targets_sim.R",
                            store = "validation")
