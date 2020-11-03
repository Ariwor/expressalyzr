peaks <- 1:8
mefl <- c(NA, 771, 2106, 6262, 15183, 45292, 136258, 291042)

mefl_dt <- data.table::data.table(Peak = peaks,
                                  MEFL = mefl
)

usethis::use_data(mefl_dt, overwrite = TRUE)
