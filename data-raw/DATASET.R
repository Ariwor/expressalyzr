peaks <- 1:8
mefl <- c(NA, 789, 1896, 4872, 15619, 47116, 143912, 333068)

mefl_dt <- data.table::data.table(Peak = peaks,
                                  MEFL = mefl
)

usethis::use_data(mefl_dt, overwrite = TRUE)
