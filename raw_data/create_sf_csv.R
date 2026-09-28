library(LBoM)
library(tidyverse)
library(readxl)
library(EMD)

year.start <- 1830

## These match for 1930 when they both exist, very good!
## For the first 3 years, no RGWR weekly birth data (data is from LBoM), so
## Olga replaced them with uniform estimation from the annual birth data
births_1881 <- 132904

births_lbom <- get_data(category = c("birth")) %>%
  filter(numdate >= 1842) %>%
  mutate(
    period_start_date = as.character(make_date(
      year.from,
      month.from,
      day.from
    )),
    period_end_date = as.character(make_date(year.to, month.to, day.to)),
    birth = birth.no.mv
  ) %>% # birth.no.mv and birth.no.uhp are equal after 1842, they use olga's estimation for 1842-1845
  select(period_start_date, period_end_date, birth) %>%
  bind_rows(data.frame(
    period_start_date = "1881-12-24",
    period_end_date = "1881-12-31",
    birth = NA
  )) %>%
  arrange(period_start_date)

births_lbom$birth[
  births_lbom$period_start_date >= "1880-12-25" &
    births_lbom$period_end_date < 1882
] <- births_1881 /
  sum(
    births_lbom$period_start_date >= "1880-12-25" &
      births_lbom$period_end_date < 1882
  ) # enter 1881 annual total births uniformly

births_post1930 <- read_excel("RGWR_London_birth_death_rows_1930-1954.xlsx") %>%
  mutate(
    period_start_date = as.character(week_ended - days(7)),
    period_end_date = as.character(week_ended),
    birth = as.numeric(live_births)
  ) %>%
  select(period_start_date, period_end_date, birth)
births <- full_join(
  births_lbom,
  births_post1930,
  by = c("period_start_date", "period_end_date"),
  relationship = "one-to-one"
) %>%
  mutate(
    births = coalesce(birth.x, birth.y)
  ) %>%
  select(period_start_date, period_end_date, births)

## Rough annual Scarlet Fever data based on LBoM package's weekly data.  Year
## starts and ends based on the first and last weeks whose dates are in that
## year, so not completely accurate.  scalet is searched as well due to a
## misspelling in one of the original data files.  Finding erroneous dates:
lbom <- read.csv("london-mort-harmonized.csv")
sf_1_strings <- lbom$cause %>% unique() %>% grep("arlet", ., value = TRUE)
sf_2_strings <- lbom$cause %>% unique() %>% grep("arlat", ., value = TRUE) # Scarlatina
## Nesting cause == "Scarlet Fever" misses a lot since they're classified as "Zymotic Diseases"
## will ignore the davenport digitized data, limited since it only has non-zero
## deaths and likely also has typos...
sf_strings <- union(sf_1_strings, sf_2_strings)

acm <- lbom %>%
  filter(
    cause == "all",
    dataset_id %in% c("acm_uk_1842-1930_age", "mort_uk_1842-1950")
  ) %>%
  arrange(period_start_date) %>%
  select(period_start_date, period_end_date, deaths) %>%
  rename(acm = deaths)

sf <- lbom %>%
  filter(
    cause %in% sf_strings,
    dataset_id != "mort_uk_1642-1845_wk_davenport"
  ) %>%
  arrange(period_start_date) %>%
  select(period_start_date, period_end_date, cause, deaths) %>%
  group_by(period_start_date, period_end_date) %>%
  summarise(deaths = sum(deaths)) %>%
  full_join(acm, by = c("period_start_date", "period_end_date")) %>%
  full_join(births, by = c("period_start_date", "period_end_date")) %>%
  arrange(period_start_date)

## add typos i found:
source("typo_fixer_helper.R")
sf_no_rgwr_typos <- fix_rgwr_typos(sf) %>%
  arrange(period_start_date) %>%
  squash_duplicates()

# Good on all fronts! Dates look continuous for rgwr period:
# anomalies <- find_period_anomalies(sf_no_rgwr_typos)
# anomalies %>% filter(bad_length | bad_connection) %>% View()
# no missing sf deaths before 1950, no missing acm:
# sf_no_rgwr_typos %>% filter(period_start_date > 1842, period_end_date < 1940) %>% pull(acm) %>% is.na() %>% which()
# no missing births after 1842:
# sf_no_rgwr_typos %>% filter(period_start_date > 1842) %>% pull(births) %>% is.na() %>% which()
# Data looks ready for use in RGWR regime!

# mutate(
#   pop = approx(x = popdat$numdate, y = popdat$pop, xout = numdate)$y
# ) %>%
# mutate(
#   birth.trend = emd(xt = birth, tt = numdate, boundary = "wave")$residue
# ) %>%
# select(numdate, date, birth, birth.trend, pop)

## population data from demographia (census + estimate from London Encyclopedia for 1939)
popdat <- data.frame(
  numdate = c(seq(1841, 1931, by = 10), decimal_date(ymd("1939-09-01"))), # Assume 1939 estimate is for sept 1, 1939, before war began for london
  pop = c(
    2207653,
    2651939,
    3188485,
    3840595,
    4713441,
    5571968,
    6506889,
    7160441,
    7386755,
    8110358,
    8615050
  )
)

post_1842_birth_emd <- sf_no_rgwr_typos %>%
  filter(period_start_date > "1842") %>%
  mutate(
    numdate = decimal_date(ymd(period_end_date)),
    birth.trend = emd(xt = births, tt = numdate, boundary = "wave")$residue
  ) %>%
  select(numdate, period_start_date, period_end_date, birth.trend)

post_1842_acm_end <- sf_no_rgwr_typos %>%
  filter(!is.na(acm)) %>%
  mutate(
    numdate = decimal_date(ymd(period_end_date)),
    acm.trend = emd(xt = acm, tt = numdate, boundary = "wave")$residue
  ) %>%
  select(numdate, period_start_date, period_end_date, acm.trend)

sf_complete <- sf_no_rgwr_typos %>%
  mutate(
    numdate = decimal_date(ymd(period_end_date)),
    pop = approx(x = popdat$numdate, y = popdat$pop, xout = numdate)$y
  ) %>%
  left_join(
    post_1842_birth_emd,
    by = c("numdate", "period_start_date", "period_end_date")
  ) %>%
  left_join(
    post_1842_acm_end,
    by = c("numdate", "period_start_date", "period_end_date")
  )

write.csv(sf_complete, "sf_complete.csv", row.names = FALSE)

## Annual Scarlet Fever data according to LBoM package.  For some reason, only
## has data with names scarlatina, scarlet.fever, and scarlet.fever.(scarlatina).
## Seems like scarlet.fever.or.scarlatina and scarlet.fever.and.streptococcal
## have weekly data but no annual data.  scarlatina and scarlet.fever columns in the annual data
## perfectly overlap between 1830 and 1841, was it to do with the change from
## LBoM to the Registrar General?  1881 is also missing annual data.
## I use the max instead of the sum here since data only ever appears in one
## column, except in the years when the columns perfectly overlap:
annual_scarlet_fever_data <- get_data(
  category = "diseaseMort",
  columns = "scarl",
  weekly = FALSE
) %>%
  dplyr::mutate(
    total.deaths = pmax(
      scarlatina,
      scarlet.fever,
      `scarlet.fever.(scarlatina)`,
      na.rm = TRUE
    )
  ) %>%
  filter(year >= year.start) # & year<=year.end)

deaths_1881 <- 2108 # looking at annual report from 1882, can see in Table 9 SF annual deaths for last 10 years, so this can be filled in!
annual_scarlet_fever_data_from_weekly <- sf_complete %>%
  mutate(year = year(ymd(period_end_date))) %>%
  group_by(year) %>%
  summarise(total.deaths = sum(deaths, na.rm = TRUE)) %>%
  filter(year >= year.start)
annual_scarlet_fever_data$total.deaths[
  annual_scarlet_fever_data$year == 1881
] <- deaths_1881
