# Web-enabled and Cross-platform PAM via Shiny

First make sure you have a very recent version of R or RStudio.

Next install required packages. Cut and paste what’s below in an R
session.

``` r
install.packages(c("pamr", "matrixStats", "superpc", "shiny", "openxlsx"))
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
```

You only need to do this once.

Then, you may run PAM any time in an R session as follows.

``` r
library(shiny)
runGitHub("PAM", "MikeJSeo")
```

That will bring up a browser window with a user interface. More details
are provided in the manual pam.pdf in this directory

Please post to the group regarding any issues. This will help us ensure
we have all the kinks ironed out before merging the code into the next
version of the `pamr` package.

Thank you,  
Michael Seo
