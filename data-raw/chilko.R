# Load packages ----------------------------------------------------------------

# install.packages("pacman")
pacman::p_load(
  dplyr,
  magrittr,
  purrr,
  readr,
  rlang,
  tibble,
  tidyr
)

# Define brood years -----------------------------------------------------------

brood_years <- c(1949:2016)

# Define stocks ----------------------------------------------------------------

stock_id <- 7
stock_name <- "Chilko"
stock_fc <- "Chilko"
stock_cu <- "Chilko"
stock_time <- "Summer & Early Summer"

stocks <- tibble::tibble(
  stock_id = stock_id,
  stock_name = stock_name,
  stock_fc = stock_fc,
  stock_cu = stock_cu,
  stock_time = stock_time
)

# Format oscillation indexes ---------------------------------------------------

oscillations <- readr::read_csv(
  paste0(
    "~/github/data-raw/sockeye/fraser/physical/2022/",
    "FC_Environmental_Data_2022_Std_Offset_NPGO_PDO.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::rename(
    brood_year = .data$BroodYear,
    pdo = .data$PDO,
    npgo_summer = .data$NPGO.Sum,
    npgo_winter = .data$NPGO.Win,
    npgo_annual = .data$NPGO.Ann
  ) %>%
  dplyr::mutate(stock_id = as.numeric(.env$stock_id)) %>%
  dplyr::relocate(.data$stock_id, .before = 1)

# Format coastal indexes -------------------------------------------------------

coastal <- readr::read_csv(
  paste0(
    "~/github/data-raw/sockeye/fraser/physical/",
    "FC_Environmental_Data_2020_Std_Offset_Real.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    brood_year = yr,
    flow = mean(c(aflow, mflow, jflow), na.rm = TRUE),
    sst_entrance = mean(c(apesst, maesst, jnesst), na.rm = TRUE),
    sst_pine = mean(c(appsst, mapsst, jnpsst, jlpsst), na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    .data$brood_year,
    .data$flow,
    .data$sst_entrance,
    .data$sst_pine
  ) %>%
  dplyr::mutate(stock_id = as.numeric(.env$stock_id)) %>%
  dplyr::relocate(.data$stock_id, .before = 1)

# Format Chilko data -----------------------------------------------------------

chilko <- read_csv(
  "~/github/data-raw/sockeye/fraser/dfo/2022/SRDATA2022.csv",
  show_col_types = FALSE
) %>%
  rename(
    stock_id = .data$PopID,
    brood_year = .data$yr,
    cycle = .data$cyc,
    r_2 = .data$rec2,
    r_3 = .data$rec3,
    r_4 = .data$rec4,
    r_5 = .data$rec5,
    recruits = .data$rec,
    spawners = .data$eff
  ) %>%
  inner_join(
    stocks,
    by = c("stock_id")
  ) %>%
  relocate(stock_name:stock_time, .before = 1) %>%
  relocate(stock_fc, .before = stock_cu) %>%
  relocate(stock_id, .before = 1) %>%
  complete(brood_year, nesting(stock_id)) %>%
  arrange(stock_id, brood_year) %>%
  select(-cycle, -juv) %>%
  ungroup() %>%
  group_by(stock_id) %>%
  mutate(
    lag_r_3 = dplyr::lag(r_3, 3),
    lag_r_4 = dplyr::lag(r_4, 4),
    lag_r_5 = dplyr::lag(r_5, 5)
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(returns = sum(lag_r_3, lag_r_4, lag_r_5, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(
    coastal,
    by = c("stock_id", "brood_year")
  ) %>%
  left_join(
    oscillations,
    by = c("stock_id", "brood_year")
  ) %>%
  dplyr::filter(.data$brood_year %in% .env$brood_years)

# Write to data/ ---------------------------------------------------------------

usethis::use_data(chilko, overwrite = TRUE)
