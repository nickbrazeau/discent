# DISCent
<!-- badges: start -->
[![master checks](https://github.com/nickbrazeau/discent/workflows/main_build/badge.svg)](https://github.com/nickbrazeau/discent/actions)
[![Codecov test coverage](https://codecov.io/gh/nickbrazeau/discent/branch/master/graph/badge.svg)](https://codecov.io/gh/nickbrazeau/discent?branch=master)
<!-- badges: end -->

<description>


## Overview
This is a simple package with a single function that runs a simple isolation by distance and identity by descent model for producing deme inbreeding spatial coefficients (DISCs). By theory, genetic relatnedness is strictly limited to measure of identity by descent. For a full mathematical overview of the _simple_ model, please see the supplements from [Brazeau et al. 2022](). 
  
Briefly, in the early 1970s, Malécot described a relationship between genetic relatedness and physical distance under the premise that pairs that are far apart are less likely to mate: isolation by distance. Capitalizing on this framework by using measures of identity by descent (IBD), we produce a deme inbreeding spatial coefficient (DISC) for point process data in continuous space. Essentially, this measure estimates the amount of "inbreeding" or within deme relatedness under the assumption that relatedness between two locations in space can be measured from the average pairwise IBD between those two locations conditional on the physical distance that seperates them. Further, we assume that geographic distance is scaled by a migration rate, which is a global parameter among all spatial locations (_i.e._ the migration rate is assumed to be constant through space). As a result, the model produces a DISC for each spatial location considered (_user defined_) and an overall migration rate parameter (_m_). Parameter optimization is achieved with a "vanilla" gradient descent algorithm. Given that "vanilla" gradient descent is used, we recommend running numerous models from a wide of array of start parameters to ensure that your final results are the global maximum and not a local maximum (_i.e._ a false peak).  


## Running the Model 
For instructions on running the model, please see our vignettes: [Vignettes](https://nickbrazeau.github.io/discent/)

## IMPORTANT NOTES
:warning: This model has planned development and is not guaranteed to be backwards compatible. Version updates that are not backwards compatible will be noted in the `NEWS.md` file.


<p align="center">
<img src="https://raw.githubusercontent.com/nickbrazeau/discent/master/R_ignore/images/discent_hex.png" width="200" height="200">
</p>
