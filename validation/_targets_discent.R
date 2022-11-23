








# run SA for each simulation type
tar_target(SAoptimparms, get_SA_wrapper_start(discdat = discdat)),

# run discent
tar_target(results, get_discentwrapper(discdat = discdat,
                                       SAoptimparms = SAoptimparms)),

# make report
tar_render(report, "validation/validation_onering.Rmd")
)




