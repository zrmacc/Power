# Add library path
library(devtools)
library(roxygen2)

# Creates package
# usethis::create_package("Power")
setwd("~/Documents/Lab/Projects/Power/Power/")

# Generates RcppExports
# pkgbuild::compile_dll(force = TRUE)

## Documents package
document()

# Install
setwd("~/Documents/Lab/Projects/Power/")
devtools::install(
  pkg = "Power", 
  reload = TRUE, 
  quick = TRUE, 
  upgrade = "never"
)

# Check package
# devtools::check()
