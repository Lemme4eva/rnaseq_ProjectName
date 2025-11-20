####For manual github setup ####
#renv::install("usethis)
#renv::install("gitcreds")
#renv::snapshot()
#library(usethis)
#library(gitcreds)
#gitcreds::gitcreds_set
#usethis::use_github



#install.packages("renv")
##renv::init()        # Initialize project environment
##renv::snapshot()    # Save package versions to renv.lock
##renv::restore()     # Recreate exact environment elsewhere


renv::install("here")
renv::install("fs")
renv::snapshot()

library(here)
library(fs)
