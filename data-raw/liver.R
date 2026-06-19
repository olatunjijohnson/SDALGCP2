# Generates data/liver.rda: a real disease-count dataset for the case-study
# tutorial. Derived from SDALGCP::PBCshp_sf (the primary biliary cirrhosis study
# of Johnson et al. 2019), trimmed to the columns a disease-mapping analysis
# needs. Re-run with:
#   Rscript data-raw/liver.R
#
# Requires the 'SDALGCP' package (the original implementation) to be installed;
# it is only needed to build the data, not to use SDALGCP2.

suppressMessages(library(sf))

stopifnot(requireNamespace("SDALGCP", quietly = TRUE))
e <- new.env()
utils::data("PBCshp_sf", package = "SDALGCP", envir = e)
pbc <- get("PBCshp_sf", envir = e)

sf::st_crs(pbc) <- 27700                      # British National Grid (OSGB eastings/northings)

liver <- pbc[, c("LSOA04CD", "pop", "IMD", "Income", "Employment", "X")]
names(liver)[names(liver) == "LSOA04CD"] <- "lsoa"
names(liver)[names(liver) == "X"]        <- "cases"
liver <- liver[, c("lsoa", "cases", "pop", "IMD", "Income", "Employment")]
rownames(liver) <- NULL

save(liver, file = "data/liver.rda", compress = "xz")
