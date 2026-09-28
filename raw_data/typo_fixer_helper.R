## this script was created with chatgpt assistance from the kz_sf_typos.csv file

fix_rgwr_typos <- function(df) {
  # Make a copy so the original data frame is not modified
  out <- df

  # Helper: update dates for a specific existing row
  fix_dates <- function(start, end, new_start = NULL, new_end = NULL) {
    idx <- out$period_start_date == start &
      out$period_end_date == end

    if (!any(idx)) {
      warning(sprintf("Could not find row: %s to %s", start, end))
      return()
    }

    if (!is.null(new_start)) {
      out$period_start_date[idx] <<- new_start
    }
    if (!is.null(new_end)) {
      out$period_end_date[idx] <<- new_end
    }
  }

  remove_dates <- function(start, end) {
    idx <- out$period_start_date == start &
      out$period_end_date == end

    if (!any(idx)) {
      warning(sprintf("Could not find row: %s to %s", start, end))
      return()
    }

    out <<- out[!idx, ]
  }

  # Helper: replace a missing row
  add_missing <- function(start, end, deaths, acm, births) {
    idx <- out$period_start_date == start &
      out$period_end_date == end

    if (!any(idx)) {
      # Row doesn't exist: add it
      out <<- rbind(
        out,
        data.frame(
          period_start_date = start,
          period_end_date = end,
          deaths = deaths,
          acm = acm,
          births = births,
          stringsAsFactors = FALSE
        )
      )

      return()
    }

    # Row exists: fill NAs, but warn about contradictions
    if (is.na(out$deaths[idx])) {
      out$deaths[idx] <<- deaths
    } else if (!is.na(deaths) && out$deaths[idx] != deaths) {
      warning(sprintf(
        "Contradiction for %s to %s: deaths is %s, correction says %s. Correction applied.",
        start,
        end,
        out$deaths[idx],
        deaths
      ))
      out$deaths[idx] <<- deaths
    }

    if (is.na(out$acm[idx])) {
      out$acm[idx] <<- acm
    } else if (!is.na(acm) && out$acm[idx] != acm) {
      warning(sprintf(
        "Contradiction for %s to %s: acm is %s, correction says %s. Correction applied.",
        start,
        end,
        out$acm[idx],
        acm
      ))
      out$acm[idx] <<- acm
    }

    if (is.na(out$births[idx])) {
      if (!is.na(births)) {
        out$births[idx] <<- births
      }
    } else if (!is.na(births) && out$births[idx] != births) {
      warning(sprintf(
        "Contradiction for %s to %s: births is %s, correction says %s. Correction applied.",
        start,
        end,
        out$births[idx],
        births
      ))
      out$births[idx] <<- births
    }
  }

  # ------------------------------------------------------------
  # Date-entry corrections
  # ------------------------------------------------------------

  fix_dates("1949-12-30", "1950-01-07", new_start = "1949-12-31")
  fix_dates("1941-11-27", "1941-12-06", new_start = "1941-11-29")
  fix_dates(
    "1941-11-20",
    "1941-11-27",
    new_start = "1941-11-22",
    new_end = "1941-11-29"
  )
  fix_dates(
    "1941-11-13",
    "1941-11-20",
    new_start = "1941-11-15",
    new_end = "1941-11-22"
  )
  fix_dates("1941-11-08", "1941-11-13", new_end = "1941-11-15")

  fix_dates("1940-07-07", "1940-07-13", new_start = "1940-07-06")
  fix_dates("1940-06-29", "1940-07-07", new_end = "1940-07-06")
  fix_dates("1939-12-29", "1940-01-06", new_start = "1939-12-30")

  fix_dates("1938-07-16", "1938-06-23", new_end = "1938-07-23")
  fix_dates("1938-06-23", "1938-07-30", new_start = "1938-07-23")

  fix_dates("1933-04-15", "1933-02-22", new_end = "1933-04-22")
  fix_dates("1933-02-22", "1933-04-29", new_start = "1933-04-22")

  fix_dates("1932-12-30", "1933-01-07", new_start = "1932-12-31")
  fix_dates("1932-08-07", "1932-08-13", new_start = "1932-08-06")
  fix_dates("1932-07-30", "1932-08-07", new_end = "1932-08-06")

  remove_dates("1927-04-23", "1927-05-30") # rm accounted for as missing
  remove_dates("1927-05-30", "1927-05-07") # rm accounted for as missing
  remove_dates("1928-04-21", "1928-05-28") # rm accounted for as missing
  remove_dates("1928-05-28", "1928-05-05") # rm accounted for as missing
  fix_dates("1928-06-23", "1928-07-30", new_end = "1928-06-30")
  fix_dates("1928-07-30", "1928-07-07", new_start = "1928-06-30")

  fix_dates("1910-01-02", "1910-01-08", new_start = "1910-01-01")

  # 1909 chain
  fix_dates("1909-04-24", "1909-05-02", new_end = "1909-05-01")
  fix_dates(
    "1909-05-02",
    "1909-05-09",
    new_start = "1909-05-01",
    new_end = "1909-05-08"
  )
  fix_dates(
    "1909-05-09",
    "1909-05-16",
    new_start = "1909-05-08",
    new_end = "1909-05-15"
  )
  fix_dates(
    "1909-05-16",
    "1909-05-23",
    new_start = "1909-05-15",
    new_end = "1909-05-22"
  )
  fix_dates(
    "1909-05-23",
    "1909-05-30",
    new_start = "1909-05-22",
    new_end = "1909-05-29"
  )
  fix_dates(
    "1909-05-30",
    "1909-06-06",
    new_start = "1909-05-29",
    new_end = "1909-06-05"
  )
  fix_dates(
    "1909-06-06",
    "1909-06-13",
    new_start = "1909-06-05",
    new_end = "1909-06-12"
  )
  fix_dates(
    "1909-06-13",
    "1909-06-20",
    new_start = "1909-06-12",
    new_end = "1909-06-19"
  )
  fix_dates(
    "1909-06-20",
    "1909-06-27",
    new_start = "1909-06-19",
    new_end = "1909-06-26"
  )
  fix_dates(
    "1909-06-27",
    "1909-07-04",
    new_start = "1909-06-26",
    new_end = "1909-07-03"
  )
  fix_dates(
    "1909-07-04",
    "1909-07-11",
    new_start = "1909-07-03",
    new_end = "1909-07-10"
  )
  fix_dates(
    "1909-07-11",
    "1909-07-18",
    new_start = "1909-07-10",
    new_end = "1909-07-17"
  )
  fix_dates(
    "1909-07-18",
    "1909-07-25",
    new_start = "1909-07-17",
    new_end = "1909-07-24"
  )
  fix_dates(
    "1909-07-25",
    "1909-08-01",
    new_start = "1909-07-24",
    new_end = "1909-07-31"
  )
  fix_dates(
    "1909-08-01",
    "1909-08-08",
    new_start = "1909-07-31",
    new_end = "1909-08-07"
  )
  fix_dates(
    "1909-08-08",
    "1909-08-15",
    new_start = "1909-08-07",
    new_end = "1909-08-14"
  )
  fix_dates(
    "1909-08-15",
    "1909-08-22",
    new_start = "1909-08-14",
    new_end = "1909-08-21"
  )
  fix_dates(
    "1909-08-22",
    "1909-08-29",
    new_start = "1909-08-21",
    new_end = "1909-08-28"
  )
  fix_dates(
    "1909-08-29",
    "1909-09-05",
    new_start = "1909-08-28",
    new_end = "1909-09-04"
  )
  fix_dates(
    "1909-09-05",
    "1909-09-12",
    new_start = "1909-09-04",
    new_end = "1909-09-11"
  )
  fix_dates(
    "1909-09-12",
    "1909-09-19",
    new_start = "1909-09-11",
    new_end = "1909-09-18"
  )
  fix_dates(
    "1909-09-19",
    "1909-09-26",
    new_start = "1909-09-18",
    new_end = "1909-09-25"
  )
  fix_dates(
    "1909-09-26",
    "1909-10-03",
    new_start = "1909-09-25",
    new_end = "1909-10-02"
  )
  fix_dates(
    "1909-10-03",
    "1909-10-10",
    new_start = "1909-10-02",
    new_end = "1909-10-09"
  )
  fix_dates(
    "1909-10-10",
    "1909-10-17",
    new_start = "1909-10-09",
    new_end = "1909-10-16"
  )
  fix_dates(
    "1909-10-17",
    "1909-10-24",
    new_start = "1909-10-16",
    new_end = "1909-10-23"
  )
  fix_dates(
    "1909-10-24",
    "1909-10-31",
    new_start = "1909-10-23",
    new_end = "1909-10-30"
  )
  fix_dates(
    "1909-10-31",
    "1909-11-07",
    new_start = "1909-10-30",
    new_end = "1909-11-06"
  )
  fix_dates(
    "1909-11-07",
    "1909-11-14",
    new_start = "1909-11-06",
    new_end = "1909-11-13"
  )
  fix_dates(
    "1909-11-14",
    "1909-11-21",
    new_start = "1909-11-13",
    new_end = "1909-11-20"
  )
  fix_dates(
    "1909-11-21",
    "1909-11-28",
    new_start = "1909-11-20",
    new_end = "1909-11-27"
  )
  fix_dates(
    "1909-11-28",
    "1909-12-05",
    new_start = "1909-11-27",
    new_end = "1909-12-04"
  )
  fix_dates(
    "1909-12-05",
    "1909-12-12",
    new_start = "1909-12-04",
    new_end = "1909-12-11"
  )
  fix_dates(
    "1909-12-12",
    "1909-12-19",
    new_start = "1909-12-11",
    new_end = "1909-12-18"
  )
  fix_dates(
    "1909-12-19",
    "1909-12-26",
    new_start = "1909-12-18",
    new_end = "1909-12-25"
  )
  fix_dates("1909-12-26", "1910-01-01", new_start = "1909-12-25")
  fix_dates(
    "1909-12-26",
    "1910-01-02",
    new_start = "1909-12-25",
    new_end = "1910-01-01"
  )

  fix_dates("1905-12-23", "1905-12-31", new_end = "1905-12-30")
  fix_dates("1905-12-31", "1906-01-06", new_start = "1905-12-30")

  fix_dates("1896-06-06", "1896-06-15", new_end = "1896-06-13")
  fix_dates(
    "1896-06-15",
    "1896-06-22",
    new_start = "1896-06-13",
    new_end = "1896-06-20"
  )
  fix_dates(
    "1896-06-22",
    "1896-06-29",
    new_start = "1896-06-20",
    new_end = "1896-06-27"
  )
  fix_dates("1896-06-29", "1896-07-04", new_start = "1896-06-27")

  fix_dates("1882-01-01", "1882-01-07", new_start = "1881-12-31")
  fix_dates("1882-12-24", "1881-12-31", new_start = "1881-12-24")
  fix_dates(
    "1882-12-24",
    "1882-01-01",
    new_start = "1881-12-24",
    new_end = "1881-12-31"
  )

  fix_dates("1874-11-21", "1874-11-27", new_end = "1874-11-28")
  fix_dates("1874-11-27", "1874-12-05", new_start = "1874-11-28")
  fix_dates("1873-11-22", "1873-11-28", new_end = "1873-11-29")
  fix_dates("1873-11-28", "1873-12-06", new_start = "1873-11-29")

  fix_dates("1868-06-06", "1868-06-18", new_end = "1868-06-13")
  fix_dates("1868-06-18", "1868-06-20", new_start = "1868-06-13")

  fix_dates("1848-10-30", "1848-10-07", new_start = "1848-09-30")
  fix_dates("1847-08-31", "1847-08-07", new_start = "1847-07-31")
  fix_dates("1847-01-23", "1847-01-31", new_end = "1847-01-30")
  fix_dates("1847-01-31", "1847-02-06", new_start = "1847-01-30")
  fix_dates("1844-04-30", "1844-04-06", new_start = "1844-03-30")
  fix_dates("1843-10-30", "1843-10-07", new_start = "1843-09-30")
  fix_dates("1842-05-30", "1842-05-07", new_start = "1842-04-30")
  fix_dates("1842-01-02", "1842-01-08", new_start = "1842-01-01")

  # ------------------------------------------------------------
  # Missing rows
  # ------------------------------------------------------------

  add_missing("1947-12-27", "1948-01-03", 0, 968, 1535)

  add_missing("1939-01-14", "1939-01-21", 0, 1274, 1054)
  add_missing("1939-01-21", "1939-01-28", 2, 1111, 969)

  add_missing("1929-11-23", "1929-11-30", 2, 1032, 1292)
  add_missing("1929-03-09", "1929-03-16", 1, 2006, 1526)
  add_missing("1929-03-16", "1929-03-23", 0, 1775, 1577)
  add_missing("1929-03-23", "1929-03-30", 2, 1331, 1347)
  add_missing("1929-03-30", "1929-04-06", 1, 1202, 1381)
  add_missing("1929-04-06", "1929-04-13", 2, 1153, 1500)
  add_missing("1929-04-13", "1929-04-20", 1, 1176, 1614)
  add_missing("1929-04-20", "1929-04-27", 0, 1092, 1391)
  add_missing("1929-04-27", "1929-05-04", 0, 1039, 1426)
  add_missing("1929-05-04", "1929-05-11", 1, 1019, 1569)

  add_missing("1928-03-10", "1928-03-17", 2, 1304, 1461)
  add_missing("1928-03-17", "1928-03-24", 0, 1361, 1581)
  add_missing("1928-03-24", "1928-03-31", 2, 1223, 1486)
  add_missing("1928-03-31", "1928-04-07", 1, 1169, 1266)
  add_missing("1928-04-07", "1928-04-14", 1, 1203, 1444)
  add_missing("1928-04-14", "1928-04-21", 2, 1137, 1524)
  add_missing("1928-04-21", "1928-04-28", 1, 1135, 1571) # acm 1135 but date typo
  add_missing("1928-04-28", "1928-05-05", 3, 1127, 1523) # acm 1127 but date typo
  add_missing("1928-05-05", "1928-05-12", 1, 948, 1554)

  add_missing("1927-03-12", "1927-03-19", 0, 1114, 1634)
  add_missing("1927-03-19", "1927-03-26", 2, 1068, 1419)
  add_missing("1927-03-26", "1927-04-02", 1, 959, 1565)
  add_missing("1927-04-02", "1927-04-09", 0, 980, 1521)
  add_missing("1927-04-09", "1927-04-16", 3, 944, 1382)
  add_missing("1927-04-16", "1927-04-23", 1, 975, 1553)
  add_missing("1927-04-23", "1927-04-30", 1, 980, 1738) # acm 980 but date typo
  add_missing("1927-04-30", "1927-05-07", 1, 959, 1526) # acm 959 but date typo
  add_missing("1927-05-07", "1927-05-14", 0, 857, 1572)

  add_missing("1908-12-26", "1909-01-02", 7, 1805, 2228)

  add_missing("1880-11-06", "1880-11-13", 84, 1636, 2538)

  add_missing("1868-12-26", "1869-01-02", 83, 1629, 2505)
  add_missing("1868-12-12", "1868-12-19", 100, 1558, 2206)
  add_missing("1867-12-14", "1867-12-21", 48, 1561, 2161)
  add_missing("1866-12-15", "1866-12-22", 30, 1377, NA)
  add_missing("1865-12-16", "1865-12-23", 43, 1590, 2109)
  add_missing("1863-12-26", "1864-01-02", 93, 1642, 2308)

  add_missing("1847-12-25", "1848-01-01", 43, 1599, 1452)

  # ------------------------------------------------------------
  # Birth, Death, ACM correction
  # ------------------------------------------------------------

  correct_row <- function(
    start,
    end,
    deaths = NULL,
    acm = NULL,
    births = NULL
  ) {
    idx <- out$period_start_date == start &
      out$period_end_date == end

    if (!any(idx)) {
      warning(sprintf("Could not find row: %s to %s", start, end))
      return()
    }

    if (!is.null(deaths)) {
      out$deaths[idx] <<- deaths
    }
    if (!is.null(acm)) {
      out$acm[idx] <<- acm
    }
    if (!is.null(births)) {
      out$births[idx] <<- births
    }
  }

  correct_row("1908-12-19", "1908-12-26", deaths = 6)
  correct_row("1868-12-19", "1868-12-26", births = 1664)
  correct_row("1890-12-27", "1891-01-03", acm = 2516)
  correct_row("1858-06-05", "1858-06-12", deaths = 47)

  out
}

squash_duplicates <- function(df) {
  df %>%
    group_by(period_start_date, period_end_date) %>%
    summarise(
      deaths = {
        x <- unique(na.omit(deaths))
        if (length(x) > 1) {
          warning(
            "Conflicting deaths values for ",
            first(period_start_date),
            " to ",
            first(period_end_date),
            ": ",
            paste(x, collapse = ", "),
            ". Keeping the first."
          )
        }
        if (length(x) == 0) NA_real_ else x[1]
      },

      acm = {
        x <- unique(na.omit(acm))
        if (length(x) > 1) {
          warning(
            "Conflicting acm values for ",
            first(period_start_date),
            " to ",
            first(period_end_date),
            ": ",
            paste(x, collapse = ", "),
            ". Keeping the first."
          )
        }
        if (length(x) == 0) NA_real_ else x[1]
      },

      births = {
        x <- unique(na.omit(births))
        if (length(x) > 1) {
          warning(
            "Conflicting births values for ",
            first(period_start_date),
            " to ",
            first(period_end_date),
            ": ",
            paste(x, collapse = ", "),
            ". Keeping the first."
          )
        }
        if (length(x) == 0) NA_real_ else x[1]
      },

      .groups = "drop"
    )
}

## data validation functions using chatgpt:
plot_date_pattern <- function(dates) {
  library(ggplot2)

  df <- data.frame(date = as.Date(dates))

  df$weekday <- factor(
    weekdays(df$date),
    levels = c(
      "Monday",
      "Tuesday",
      "Wednesday",
      "Thursday",
      "Friday",
      "Saturday",
      "Sunday"
    )
  )

  ggplot(df, aes(x = date, y = weekday, colour = weekday)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
      x = "Date",
      y = NULL,
      colour = "Weekday"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

find_period_anomalies <- function(df) {
  start <- as.Date(df$period_start_date)
  end <- as.Date(df$period_end_date)

  # 1. Period is not exactly 7 days
  bad_length <- (end - start) != 7

  # 2. End of row k != start of row k+1
  bad_connection <- end[-nrow(df)] != start[-1]

  # 3. Neither start nor end is repeated on the next row
  bad_nonrepeat <- (end[-nrow(df)] != end[-1] &
    start[-nrow(df)] != start[-1])

  # Pad conditions 2 and 3 so they align with rows
  bad_connection <- c(bad_connection, FALSE)
  bad_nonrepeat <- c(bad_nonrepeat, FALSE)

  # Return indices, plus which rule(s) were violated
  anomalies <- which(
    bad_length | bad_connection | bad_nonrepeat
  )

  data.frame(
    index = anomalies,
    bad_length = bad_length[anomalies],
    bad_connection = bad_connection[anomalies],
    bad_nonrepeat = bad_nonrepeat[anomalies]
  )
}
