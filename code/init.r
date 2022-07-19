#' Todo
#' 
#' - read covid 8, 9
#' - merge 8, 9
#' - merge 1:7
#' 
#' # Visits
#' 
#' - Binomial model
#' - Covariate before covid -> predict prob drop out during covid
#' 
#' # Clustering - attitude to vaccine

library(tidyverse)
library(haven)
here::i_am("code/init.r")

devtools::load_all('~/Code/R/ktools/')
read_in_zip <- function(file, pattern, readfn = haven::read_dta) {
    name <- find_in_zip(file, pattern)
    unar(file, exdir = tempdir(), force_directory = F)
    readfn(paste0(tempdir(), '/', name))
}

#' # Wave 9
#'
file <- "~/Seafile/SHARE/sharew9ca_rel8-0-0_ALL_datasets_stata.zip"
list_zip(file)

# coverscreen: meta data
cvr <- read_in_zip(file, "_cv_r")
names(cvr)
cvr %>% count(interview)
cvr %>% count(waveid, firstwave)
cvr %>% count(interview, waveid) %>% print(n = Inf)

# end of life
xt <- haven::read_dta(file)
#' there is covid death
xt %>% count(xt011_)

# covid
covid <- read_in_zip(file, "-0_ca")
covid %>% names

which(xt$mergeid %in% cvr$mergeid) %>% length()

list_zip(file)
