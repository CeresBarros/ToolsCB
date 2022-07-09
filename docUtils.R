library(Require)

Require("roxygen2md")

## run first to set up markdown support;
## then stage changes in git to see subsequent changes
roxygen2md(scope = "none")

## deals only with \code{blah} to `blah` changes
roxygen2md(scope = "simple")

## deals with most other changes
roxygen2md(scope = "full")

## show remaining changes that will need to be done manually
find_rd()
