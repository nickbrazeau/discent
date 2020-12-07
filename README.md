# DISCent
<!-- badges: start -->
[![Travis build status](https://travis-ci.com/nickbrazeau/discent.svg?branch=main)](https://travis-ci.com/nickbrazeau/discent)
<!-- badges: end -->

<description>

## IMPORTANT NOTES
:warning: This model has planned development and is guaranteed to be backwards compatible. Version updates that are not backwards compatible will be noted in the `NEWS.md` file.


## Overview
This is a simple package with a single function that runs a simple isolation by distance and identity by descent model for producing deme inbreeding spatial coefficients (DISCs). By theory, genetic relatnedness is strictly limited to measure of identity by descent. For a full mathematical overview of the _simple_ model, please see [TODO](). 
  
Briefly, in the early 1970s, Malécot described a relationship between genetic relatedness and physical distance under the premise that pairs that are far apart are less likely to mate: isolation by distance. Capitalizing on this framework by using measures of identity by descent (IBD), we produce a deme inbreeding spatial coefficient (DISC) for point process data in continuous space. Essentially, this measure estimates the amount of "inbreeding" or within deme relatedness under the assumption that relatedness between two locations in space can be measured from the average pairwise IBD between those two locations conditional on the physical distance that seperates them. Further, we assume that geographic distance is scaled by a migration rate, which is a global parameter among all spatial locations (_i.e._ the migration rate is assumed to be constant through space). As a result, the model produces a DISC for each spatial location considered (_user defined_) and an overall migration rate parameter (_m_). Parameter optimization is achieved with a "vanilla" gradient descent algorithm. Given that "vanilla" gradient descent is used, we recommend running numerous models from a wide of array of start parameters to ensure that your final results are the global maximum and not a local maximum (_i.e._ a false peak).  



## Installation 
```
install.packages("remotes")
remotes::install_github("nickbrazeau/discent")
```
### Dependencies
`discent` relies on the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package, which requires certain OS-specific dependencies: 

* Windows
    - Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
* Mac OS X
    - Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    - Within XCode go to Preferences : Downloads and install the Command Line Tools
* Linux (Debian/Ubuntu)
    - Install the core software development utilities required for R package development as well as LaTeX by executing
    ```
    sudo apt-get install r-base-dev texlive-full
    ```

Assuming everything installed correctly, you can now load the package:

```
library(discent)
```


## Running the Model 
For instructions on running the model, please see the vignette: [Running the Model]()
