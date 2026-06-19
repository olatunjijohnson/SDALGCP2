# Generates data/liver.rda from the frozen source in data-raw/liver_source.rds.
# Re-run with:  Rscript data-raw/liver.R
#
# liver_source.rds is a one-time, self-contained snapshot of the study data
# (incident primary biliary cirrhosis counts and area deprivation covariates by
# LSOA, North East England; Johnson, Diggle and Giorgi 2019); regenerating the
# dataset does not require any other package.

suppressMessages(library(sf))

src <- readRDS("data-raw/liver_source.rds")     # sf, EPSG:27700

liver <- src
names(liver)[names(liver) == "LSOA04CD"] <- "lsoa"
names(liver)[names(liver) == "X"]        <- "cases"
liver <- liver[, c("lsoa", "cases", "pop", "IMD", "Income", "Employment")]
rownames(liver) <- NULL

save(liver, file = "data/liver.rda", compress = "xz")
