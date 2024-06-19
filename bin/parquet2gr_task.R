
library(devtools)
library(data.table)
devtools::load_all('~/git/chromunity')
args = commandArgs(trailingOnly = TRUE)
concatemers = parquet2gr(args[1])
saveRDS(concatemers, 'all_concatemers.rds')
