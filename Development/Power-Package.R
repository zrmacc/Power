# Add library path
library(devtools)
library(roxygen2)

# Creates package
base_dir <- "~/Documents/Lab/Projects/Power/"
pkg_dir <- file.path(base_dir, "Power")

# Generates RcppExports
# pkgbuild::compile_dll(force = TRUE)

# Document package
setwd(pkg_dir)
document()

# Install
setwd(base_dir)
devtools::install(pkg = "Power", reload = TRUE)

# Test.
setwd(pkg_dir)
devtools::test()

# Check package
# devtools::check()
