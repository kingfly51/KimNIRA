# Shared time-window construction. Never compress missing observations before
# constructing lags. The audit contains original input row numbers.
.kim_time_pairs <- function(data, vars, id, time = NULL, lags = 1L, min_obs = 1L,
                            max_time_gap = NULL, min_time_gap = NULL,
                            consec = NULL, dayvar = NULL, beepvar = NULL,
                            allow_overnight = FALSE, time_unit = "auto",
                            beeps_per_day = NULL) {
  # Day/beep is an alternative to a global index, not extra required metadata.
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  day_beep <- !is.null(dayvar) && !is.null(beepvar)
  if (is.null(time) && !day_beep) {
    time <- if ("time" %in% names(data)) "time" else if (!is.null(consec)) consec else NULL
    if (is.null(time)) stop("Supply time (original measurement index) OR dayvar + beepvar.")
  }
  columns <- c(vars, id, time, consec, dayvar, beepvar)
  if (!is.data.frame(data) || !all(columns %in% names(data)))
    stop("data must contain all specified node, ID and time columns.")
  if (!length(vars) || anyDuplicated(vars) ||
      !all(vapply(data[vars], function(x) is.numeric(x) && all(is.na(x) | x %in% 0:1), logical(1))))
    stop("Node columns must be distinct numeric binary (0/1) variables.")
  if (anyNA(data[c(id, time, consec, dayvar, beepvar)]))
    stop("Missing IDs or time/schedule metadata are not allowed.")
  if (!is.null(time) && anyDuplicated(data[c(id, time)])) stop("Repeated person-time indices.")
  if (day_beep && anyDuplicated(data[c(id, dayvar, beepvar)])) stop("Repeated person-day-beep indices.")
  if (!length(lags) || any(!is.finite(lags) | lags < 1 | lags != floor(lags)))
    stop("lags must be positive integers.")
  if (length(allow_overnight) != 1L || is.na(allow_overnight) || !is.logical(allow_overnight))
    stop("allow_overnight must be TRUE or FALSE.")
  if (!is.null(beepvar) && is.null(dayvar)) stop("beepvar requires dayvar.")
  if (!is.null(beeps_per_day) && (length(beeps_per_day) != 1L ||
      !is.numeric(beeps_per_day) || !is.finite(beeps_per_day) ||
      beeps_per_day < 1 || beeps_per_day != floor(beeps_per_day)))
    stop("beeps_per_day must be a positive integer from the planned schedule.")
  if (!is.null(beeps_per_day) && !day_beep) stop("beeps_per_day requires dayvar + beepvar.")
  for (nm in c(consec, beepvar)) {
    x <- data[[nm]]
    if (!is.numeric(x) || any(!is.finite(x) | x != floor(x)))
      stop("consec and beepvar columns must contain finite integer schedule indices.")
  }
  for (bound in list(min_time_gap, max_time_gap))
    if (!is.null(bound) && (length(bound) != 1L || !is.numeric(bound) ||
        !is.finite(bound) || bound < 0)) stop("Time-gap bounds must be finite nonnegative numbers.")
  if (!is.null(min_time_gap) && !is.null(max_time_gap) && min_time_gap > max_time_gap)
    stop("min_time_gap cannot exceed max_time_gap.")
  time_unit <- match.arg(time_unit, c("auto", "secs", "mins", "hours", "days", "weeks"))
  tt <- if (is.null(time)) NULL else data[[time]]
  dated <- inherits(tt, "POSIXt") || inherits(tt, "Date")
  if (!is.null(tt) && !dated && !is.numeric(tt)) stop("time must be numeric, Date or POSIXct/POSIXlt.")
  if (is.null(tt) && (!is.null(min_time_gap) || !is.null(max_time_gap) || time_unit != "auto"))
    stop("Elapsed-time bounds require an explicit time column; day/beep alone records the schedule, not elapsed time.")
  if (!dated && time_unit != "auto") stop("Numeric time uses its own units; use time_unit = 'auto'.")
  unit <- if (time_unit == "auto") {
    if (is.null(tt)) "unavailable" else if (inherits(tt, "Date")) "days" else if (dated) "secs" else "numeric"
  } else time_unit
  if (any(!is.finite(as.numeric(tt)))) stop("Non-finite time values.")
  schedule <- if (!is.null(consec)) data[[consec]] else NULL
  # A numeric integer time column means ORIGINAL scheduled indices. With a
  # day/beep specification it can also supply the global overnight index.
  if (is.null(schedule) && !is.null(tt) && !dated && all(tt == floor(tt)))
    schedule <- tt
  if (day_beep) {
    days <- data[[dayvar]]
    if (!inherits(days, "Date") && (!is.numeric(days) || any(!is.finite(days) | days != floor(days))))
      stop("dayvar must be an integer study-day index or Date, preserving missing days.")
    if (any(data[[beepvar]] < 1)) stop("beepvar must use positive scheduled indices starting at one.")
    if (!is.null(beeps_per_day) && any(data[[beepvar]] > beeps_per_day))
      stop("Observed beep exceeds beeps_per_day.")
    if (allow_overnight && is.null(schedule) && is.null(beeps_per_day))
      stop("For overnight day/beep transitions, supply beeps_per_day (planned count) or a global integer time index. The count is never inferred from observed beeps.")
    if (allow_overnight && is.null(schedule)) {
      schedule <- as.numeric(days) * beeps_per_day + data[[beepvar]] - 1
    }
  }
  if (is.null(schedule) && !day_beep && is.null(max_time_gap))
    warning("No schedule continuity or maximum time gap supplied: adjacent recorded rows are assumed consecutive. Supply consec, dayvar/beepvar, or max_time_gap.", call. = FALSE)
  depth <- max(lags)
  audits <- list(); windows <- list(); a <- 0L
  for (pid in unique(as.character(data[[id]]))) {
    rows <- which(as.character(data[[id]]) == pid)
    rows <- if (is.null(tt)) rows[order(data[[dayvar]][rows], data[[beepvar]][rows])] else rows[order(tt[rows])]
    n <- length(rows)
    if (n <= depth) next
    d <- data[rows, , drop = FALSE]
    gap <- if (is.null(tt)) rep(NA_real_, n - 1L) else if (dated)
      as.numeric(difftime(d[[time]][-1L], d[[time]][-n], units = unit)) else diff(d[[time]])
    bad_schedule <- rep(FALSE, n - 1L)
    cross_day <- if (is.null(dayvar)) rep(FALSE, n - 1L) else d[[dayvar]][-1L] != d[[dayvar]][-n]
    if (!is.null(schedule)) {
      dc <- diff(schedule[rows])
      if (any(dc <= 0)) stop("Global measurement indices must increase strictly within each person in time order.")
      bad_schedule <- dc != 1
    }
    if (!is.null(beepvar)) {
      db <- diff(d[[beepvar]])
      if (any(!cross_day & db <= 0)) stop("beepvar must increase strictly within each person-day.")
      bad_schedule <- bad_schedule | (!cross_day & db != 1)
    }
    bad_gap <- if (is.null(tt)) rep(FALSE, n - 1L) else gap <= 0
    if (!is.null(min_time_gap)) bad_gap <- bad_gap | gap < min_time_gap
    if (!is.null(max_time_gap)) bad_gap <- bad_gap | gap > max_time_gap
    end <- seq.int(depth + 1L, n)
    ix <- vapply(0:depth, function(l) rows[end - l], integer(length(end)))
    dim(ix) <- c(length(end), depth + 1L)
    reason <- rep("retained", length(end))
    for (j in seq_along(end)) {
      edges <- seq.int(end[j] - depth, end[j] - 1L)
      used <- ix[j, c(1L, lags + 1L)]
      reason[j] <- if (any(bad_schedule[edges])) "nonconsecutive" else
        if (!allow_overnight && any(cross_day[edges])) "overnight" else
        if (any(bad_gap[edges])) "time_gap" else
        if (anyNA(data[used, vars, drop = FALSE])) "missing" else "retained"
    }
    if (sum(reason == "retained") < min_obs) reason[reason == "retained"] <- "min_obs"
    a <- a + 1L
    audits[[a]] <- data.frame(person = pid, position = end, response_row = ix[, 1L],
      first_row = ix[, depth + 1L], elapsed = vapply(seq_along(end), function(j)
        sum(gap[seq.int(end[j] - depth, end[j] - 1L)]), numeric(1)), reason = reason,
      stringsAsFactors = FALSE)
    windows[[a]] <- ix
  }
  if (!length(audits)) stop("No candidate lag windows; joint lag-1/lag-2 requires three ordered observations.")
  audit <- do.call(rbind, audits); ix <- do.call(rbind, windows)
  keep <- audit$reason == "retained"
  levels <- c("retained", "nonconsecutive", "overnight", "time_gap", "missing", "min_obs")
  counts <- table(factor(audit$reason, levels = levels))
  list(X = lapply(lags, function(l) as.matrix(data[ix[keep, l + 1L], vars, drop = FALSE])),
       Y = as.matrix(data[ix[keep, 1L], vars, drop = FALSE]),
       person = audit$person[keep], position = audit$position[keep],
       rows = ix[keep, , drop = FALSE], audit = audit,
       summary = data.frame(reason = levels, n = as.integer(counts)), time_unit = unit,
       excluded = unique(audit$person[audit$reason == "min_obs"]))
}
