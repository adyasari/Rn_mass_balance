# **************************
# project setup ----
# author: peter fuleky
# **************************

# remove all objects from global environment
rm(list = ls())
# rm(list = ls(pattern = glob2rx("*__*")))

# detach all loaded packages
if (!is.null(names(sessionInfo()$otherPkgs))) {
invisible(
  suppressMessages(
    suppressWarnings(
      lapply(
        paste("package:", names(sessionInfo()$otherPkgs), sep=""), 
        detach, 
        character.only = TRUE, 
        unload = TRUE
        )
      )
    )
  )
}

# set location relative to project root. use here() instead of path names
here::i_am("R/setup.R")

# load necessary packages
library("here")
library("conflicted")
library("tidyverse")

# detect conflicts across packages and assign preferences
conflict_scout()
conflict_prefer("filter", "dplyr") # dplyr v stats
conflict_prefer("lag", "dplyr") # dplyr v stats

# top level project directory
here()

# store sensitive info in project specific .Renviron file (must end with \n)
# Sys.getenv("api_username")
# Sys.getenv("api_key")

# load user defined utility functions
source(here("R", "util_funs.R"))

