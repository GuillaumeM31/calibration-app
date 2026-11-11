# =============================================================================
# Calibration Curve Analysis - v3.5
# Weighted regression, quantification, diagnostics and group statistics
#
# Main changes vs v3.4 (details in the "Guide" tab):
#  - Bug fixes: manual table editing (0-based DT columns), Cook's distance,
#    blank (C = 0) handling in accuracy and weighting, stale results after
#    re-fitting, manual sample entry never reaching the statistics tab,
#    3-SD filter creating empty rows, equation formatting.
#  - "Intelligent" point selection now brackets the sample range and never
#    removes standards silently (influential points are flagged only).
#  - Everything is reactive: no "Calculate" buttons, nothing can go stale.
#  - Responsive layout: plots resize with the window, larger fonts,
#    per-tab collapsible sidebars, full-screen cards.
#  - New: LOD/LOQ, 95% CI on concentrations (OLS linear), Mandel test,
#    saturation warning, below/above range status, concentration unit,
#    dilution factor, non-parametric tests, robust outlier flags.
# =============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
library(dplyr)

# =============================================================================
# Helpers: file reading and formatting
# =============================================================================

# Detect the separator (; , or tab) and read a CSV file
detect_and_read_csv <- function(file_path) {
  first_line <- readLines(file_path, n = 1, warn = FALSE, encoding = "UTF-8")
  count <- function(pattern) lengths(regmatches(first_line, gregexpr(pattern, first_line)))
  separators <- c(";" = count(";"), "," = count(","), "\t" = count("\t"))
  detected_sep <- if (max(separators) == 0) ";" else names(which.max(separators))

  read_with <- function(sep) {
    read.csv(file_path, sep = sep, stringsAsFactors = FALSE,
             fileEncoding = "UTF-8-BOM", check.names = FALSE)
  }
  df <- read_with(detected_sep)
  if (ncol(df) == 1) {
    for (sep in setdiff(c(";", ",", "\t"), detected_sep)) {
      df_test <- tryCatch(read_with(sep), error = function(e) NULL)
      if (!is.null(df_test) && ncol(df_test) > 1) {
        df <- df_test
        detected_sep <- sep
        break
      }
    }
  }
  names(df) <- trimws(names(df))
  attr(df, "separator") <- switch(detected_sep, ";" = "semicolon (;)",
                                  "," = "comma (,)", "\t" = "tab")
  df
}

# Convert text to numbers, accepting decimal commas and spaces as thousand separators
to_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x <- gsub("[[:space:] ]", "", x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

# Format numbers with significant digits (works for 1e-6 as well as 1e8)
fmt <- function(x, digits = 4) {
  out <- formatC(signif(x, digits), digits = digits, format = "g")
  out[!is.finite(x)] <- "NA"
  trimws(out)
}

# "+ 12.3" / "- 12.3" for equation terms
signed <- function(x, digits = 4) {
  paste(if (is.finite(x) && x < 0) "−" else "+", fmt(abs(x), digits))
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) "< 0.0001" else sprintf("%.4f", p)
}

p_stars <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

# Common ggplot theme (large fonts, readable on projector and laptop)
theme_app <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.title.position = "plot",
      plot.subtitle = element_text(colour = "#555555"),
      plot.caption = element_text(colour = "#666666", size = rel(0.75)),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.margin = margin(10, 15, 10, 10)
    )
}

# Concentration axis: pseudo-log scale when the range spans > 2 decades.
# (pseudo-log keeps the blank at 0 visible, unlike a pure log scale)
conc_x_scale <- function(x, name = waiver()) {
  pos <- x[is.finite(x) & x > 0]
  if (length(pos) < 2 || max(pos) / min(pos) < 100) {
    return(scale_x_continuous(name = name))
  }
  brks <- scales::breaks_log(6)(range(pos))
  if (any(x == 0, na.rm = TRUE)) brks <- c(0, brks)
  scale_x_continuous(
    name = name,
    trans = scales::pseudo_log_trans(sigma = min(pos) / 2, base = 10),
    breaks = brks,
    labels = function(b) fmt(b, 3)
  )
}

# =============================================================================
# Helpers: calibration calculations
# =============================================================================

# Weights normalised to mean 1 (does not change the fit, keeps sigma interpretable).
# The blank (C = 0) gets the weight of the lowest standard instead of an
# arbitrary huge weight (v3.4 used 10x the maximum, forcing the line through 0).
compute_weights <- function(conc, scheme) {
  if (scheme == "none") return(NULL)
  pos <- conc[conc > 0]
  if (!length(pos)) return(NULL)
  x <- ifelse(conc > 0, conc, min(pos))
  w <- if (scheme == "1/x") 1 / x else 1 / x^2
  w / mean(w)
}

# Intelligent selection: keep the calibration levels whose responses cover the
# sample responses, plus ONE standard below and ONE above (bracketing), so that
# every sample is interpolated, never extrapolated.
select_intelligent <- function(df, sample_resp, margin, include_zero, min_points) {
  n <- nrow(df)
  sample_resp <- sample_resp[is.finite(sample_resp) & sample_resp > 0]
  if (!length(sample_resp)) {
    return(list(used = rep(TRUE, n),
                note = paste("No sample data loaded yet: all calibration points are used.",
                             "Load samples (tab 2) to activate intelligent selection.")))
  }
  std <- which(df$Concentration > 0)
  std <- std[order(df$Concentration[std])]
  if (length(std) < 2) {
    return(list(used = rep(TRUE, n), note = "Fewer than 2 non-zero standards: all points used."))
  }
  r <- df$Response[std]
  lo <- min(sample_resp) * (1 - margin)
  hi <- max(sample_resp) * (1 + margin)

  inside <- which(r >= lo & r <= hi)
  below <- which(r < lo)
  above <- which(r > hi)
  candidates <- c(inside,
                  if (length(below)) max(below),
                  if (length(above)) min(above))
  k <- seq(min(candidates), max(candidates))

  n_zero <- if (include_zero) sum(df$Concentration == 0) else 0
  while (length(k) + n_zero < min_points && (min(k) > 1 || max(k) < length(std))) {
    if (min(k) > 1) k <- c(min(k) - 1, k)
    if (length(k) + n_zero < min_points && max(k) < length(std)) k <- c(k, max(k) + 1)
  }

  used <- rep(FALSE, n)
  used[std[k]] <- TRUE
  if (include_zero) used[df$Concentration == 0] <- TRUE

  note <- sprintf(
    paste("Sample responses: %s to %s (window ±%d%%: %s to %s).",
          "Selected levels: %s to %s, including one bracketing standard on each side where available.%s"),
    fmt(min(sample_resp), 3), fmt(max(sample_resp), 3), round(margin * 100),
    fmt(lo, 3), fmt(hi, 3),
    fmt(df$Concentration[std[min(k)]], 4), fmt(df$Concentration[std[max(k)]], 4),
    if (!length(above)) " Warning: no standard above the highest sample." else ""
  )
  list(used = used, note = note)
}

# Fit one calibration model. Returns a list used everywhere else.
fit_calibration <- function(d, type, scheme) {
  notes <- character(0)
  d <- d[, c("Concentration", "Response")]
  if (type %in% c("log", "power")) {
    keep <- d$Concentration > 0 & (type != "power" | d$Response > 0)
    if (any(!keep)) {
      notes <- c(notes, sprintf("%d point(s) with zero concentration/response dropped (not allowed in the %s model).",
                                sum(!keep), type))
    }
    d <- d[keep, , drop = FALSE]
  }
  n_par <- if (type == "quadratic") 3 else 2
  if (nrow(d) < n_par + 1) {
    stop(sprintf("At least %d calibration points are needed for a %s model (%d available).",
                 n_par + 1, type, nrow(d)))
  }
  w <- compute_weights(d$Concentration, scheme)
  if (type == "power" && !is.null(w)) {
    w <- NULL
    notes <- c(notes, "Weighting is ignored for the power model: the log-log fit already gives relative (proportional) errors.")
  }
  d$.w <- if (is.null(w)) 1 else w

  fml <- switch(type,
    linear    = Response ~ Concentration,
    log       = Response ~ log(Concentration),
    quadratic = Response ~ Concentration + I(Concentration^2),
    power     = log(Response) ~ log(Concentration)
  )
  model <- lm(fml, data = d, weights = .w)

  list(model = model, type = type,
       scheme = if (is.null(w)) "none" else scheme,
       weighted = !is.null(w), data = d, notes = notes)
}

predict_response <- function(f, conc) {
  p <- suppressWarnings(predict(f$model, newdata = data.frame(Concentration = conc)))
  if (f$type == "power") exp(p) else unname(p)
}

# Inverse prediction: response -> concentration
conc_from_response <- function(f, y) {
  b <- unname(coef(f$model))
  x <- switch(f$type,
    linear = (y - b[1]) / b[2],
    log    = exp((y - b[1]) / b[2]),
    power  = suppressWarnings(exp((log(y) - b[1]) / b[2])),
    quadratic = {
      cc <- b[1] - y
      D <- b[2]^2 - 4 * b[3] * cc
      out <- rep(NA_real_, length(y))
      ok <- is.finite(D) & D >= 0
      # Numerically stable root on the increasing branch (also valid when b[3] ~ 0)
      out[ok] <- -2 * cc[ok] / (b[2] + sqrt(D[ok]))
      out
    }
  )
  x[!is.finite(x)] <- NA_real_
  unname(x)
}

# 95% confidence interval half-width on an interpolated concentration
# (classical formula, single measurement; OLS linear model only)
ci_halfwidth <- function(f, y0, level = 0.95) {
  if (f$type != "linear" || f$weighted) return(NULL)
  x <- f$data$Concentration
  y <- f$data$Response
  n <- length(x)
  b <- unname(coef(f$model)[2])
  s <- summary(f$model)$sigma
  sx0 <- (s / abs(b)) * sqrt(1 + 1 / n + (y0 - mean(y))^2 / (b^2 * sum((x - mean(x))^2)))
  qt(1 - (1 - level) / 2, n - 2) * sx0
}

# LOD / LOQ from the residual SD at the low end of the curve (ICH Q2 approach).
# For weighted fits, SD at the lowest standard = sigma_w / sqrt(w_lowest).
lod_loq <- function(f) {
  if (!f$type %in% c("linear", "quadratic")) return(NULL)
  slope <- unname(coef(f$model)[2])
  pos <- which(f$data$Concentration > 0)
  if (!length(pos) || !is.finite(slope) || slope <= 0) return(NULL)
  lowest <- pos[which.min(f$data$Concentration[pos])]
  sd_low <- summary(f$model)$sigma / sqrt(f$data$.w[lowest])
  list(lod = 3.3 * sd_low / slope, loq = 10 * sd_low / slope)
}

# Mandel-type test: is the quadratic term significant? (linear models only)
mandel_p <- function(f) {
  if (f$type != "linear" || nrow(f$data) < 4) return(NA_real_)
  q <- lm(Response ~ Concentration + I(Concentration^2), data = f$data, weights = .w)
  tryCatch(anova(f$model, q)$`Pr(>F)`[2], error = function(e) NA_real_)
}

# Apply DT cell edits. DT returns 0-based column indices when rownames = FALSE
# (v3.4 used them as 1-based, so edits went to the wrong column).
apply_edits <- function(df, info, numeric_cols) {
  for (k in seq_len(nrow(info))) {
    i <- info$row[k]
    j <- info$col[k] + 1
    if (j < 1 || j > ncol(df)) next
    col <- names(df)[j]
    df[i, j] <- if (col %in% numeric_cols) to_numeric(info$value[k]) else as.character(info$value[k])
  }
  df
}

add_empty_row <- function(df) {
  df[nrow(df) + 1, ] <- NA
  chr <- vapply(df, is.character, logical(1))
  if (any(chr)) df[nrow(df), chr] <- ""
  df
}

# PDF device: cairo (better Unicode) on desktop R, base pdf() in the browser (Shinylive/webR)
pdf_device <- function() if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"

STATUS_LEVELS <- c("✓ In range", "⚠ Below range", "⚠ Above range", "✗ Not calculable")
STATUS_COLOURS <- c("#d4edda", "#fff3cd", "#f8d7da", "#e2e3e5")
SAMPLE_NUMERIC <- c("Response", "Theoretical_Concentration", "Dilution_Factor")

# =============================================================================
# UI
# =============================================================================

app_css <- "
  .dataTables_wrapper, table.dataTable { font-size: 0.95rem; }
  .card-header { font-weight: 600; font-size: 1.05rem; }
  .bslib-value-box .value-box-title { font-size: 0.95rem; margin-bottom: 0.2rem; }
  .bslib-value-box .value-box-value { font-size: 1.3rem; }
  .bslib-value-box .value-box-area { padding: 0.6rem 1rem !important; }
  #shiny-notification-panel { width: 380px; top: 70px; bottom: unset; }
  .metric-table td:first-child { color: #555; width: 55%; }
  .metric-table td { padding: 0.35rem 0.5rem; }
  .eq { font-size: 1.15rem; font-weight: 600; color: #1e7e46; font-family: 'Courier New', monospace; }
  .note { font-size: 0.9rem; color: #555; }
  .sidebar-title { font-weight: 600; }
"

plot_card <- function(title, output_id, height = "48vh", min_height = "360px",
                      header_extra = NULL, footer = NULL) {
  card(
    full_screen = TRUE, height = height, min_height = min_height,
    card_header(class = "d-flex justify-content-between align-items-center",
                title, header_extra),
    card_body(plotOutput(output_id), padding = 6),
    if (!is.null(footer)) card_footer(footer)
  )
}

two_cols <- function(...) {
  layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)), ...)
}

ui <- page_navbar(
  title = "Calibration & Quantification v3.5",
  window_title = "Calibration v3.5",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50",
                   font_scale = 1.1),
  fillable = FALSE,
  header = tags$head(tags$style(HTML(app_css))),

  # ---------------------------------------------------------------- Tab 1
  nav_panel(
    "1. Calibration",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340, title = "Calibration setup",
        fileInput("calib_file", "Upload calibration CSV", accept = c(".csv", ".txt")),
        uiOutput("example_ui"),
        helpText("Two columns: concentration and response (area/signal).",
                 "Separators ; , or tab and decimal commas are detected automatically.",
                 "You can also type or correct values in the table below the plot."),
        textInput("conc_unit", "Concentration unit", value = "",
                  placeholder = "e.g. nM, µg/mL"),
        selectInput("regression_type", "Regression model",
                    choices = c("Linear" = "linear", "Quadratic" = "quadratic",
                                "Logarithmic" = "log", "Power (log-log)" = "power")),
        selectInput("weighting_scheme", "Weighting",
                    choices = c("None (OLS)" = "none", "1/X" = "1/x", "1/X²" = "1/x2")),
        helpText("Use 1/X or 1/X² when the error grows with concentration (typical for LC-MS",
                 "and wide ranges). The blank gets the weight of the lowest standard."),
        radioButtons("point_selection_mode", "Calibration points",
                     choices = c("Use all points" = "all",
                                 "Bracket the samples (intelligent)" = "intelligent",
                                 "Manual selection" = "manual"),
                     selected = "intelligent"),
        conditionalPanel(
          "input.point_selection_mode == 'intelligent'",
          sliderInput("safety_margin", "Margin around sample responses (%)",
                      min = 0, max = 50, value = 10, step = 5),
          checkboxInput("include_zero", "Include blank (C = 0)", TRUE),
          sliderInput("min_points", "Minimum number of points",
                      min = 3, max = 10, value = 5, step = 1)
        ),
        conditionalPanel(
          "input.point_selection_mode == 'manual'",
          uiOutput("manual_point_selector")
        ),
        checkboxInput("drop_influential",
                      "Automatically exclude influential points (Cook's D > 4/n) - not recommended",
                      FALSE)
      ),
      uiOutput("calib_warnings"),
      layout_column_wrap(
        width = "200px", fill = FALSE,
        value_box(title = textOutput("vb_r2_label", inline = TRUE),
                  value = textOutput("vb_r2"), theme = "primary"),
        value_box(title = "Points used", value = textOutput("vb_points"), theme = "secondary"),
        value_box(title = "Working range", value = textOutput("vb_range"), theme = "secondary"),
        value_box(title = "LOD / LOQ", value = textOutput("vb_lod"), theme = "secondary")
      ),
      plot_card(
        "Calibration curve", "calib_plot", height = "62vh", min_height = "420px",
        header_extra = downloadButton("download_calib_plot", "PDF",
                                      class = "btn-sm btn-outline-primary"),
        footer = div(
          class = "d-flex flex-wrap gap-4",
          checkboxInput("show_excluded_points", "Show unused points", FALSE),
          checkboxInput("show_samples_on_curve", "Show samples on the curve", TRUE),
          checkboxInput("calib_log_axes", "Log-log axes", FALSE)
        )
      ),
      two_cols(
        card(card_header("Equation and point selection"),
             uiOutput("regression_equation"),
             uiOutput("point_selection_summary")),
        card(card_header("Calibration data (editable: double-click a cell)"),
             div(class = "d-flex gap-2 mb-2",
                 actionButton("add_calib_row", "Add row", class = "btn-sm btn-success"),
                 actionButton("remove_calib_row", "Remove last row", class = "btn-sm btn-outline-danger")),
             DTOutput("calib_table"))
      )
    )
  ),

  # ---------------------------------------------------------------- Tab 2
  nav_panel(
    "2. Samples",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340, title = "Sample data",
        fileInput("sample_file", "Upload sample CSV", accept = c(".csv", ".txt")),
        helpText("Required: sample name + response (area).",
                 "Optional: Theoretical_Concentration (accuracy), Dilution_Factor,",
                 "and any text column (treatment, class...) used for grouping."),
        div(class = "d-flex gap-2",
            actionButton("add_sample_row", "Add row", class = "btn-sm btn-success"),
            actionButton("remove_sample_row", "Remove last row", class = "btn-sm btn-outline-danger")),
        hr(),
        radioButtons("csv_format", "Export format",
                     choices = c("Excel FR (; and decimal comma)" = "fr",
                                 "International (, and decimal point)" = "en")),
        downloadButton("download_results", "Download results", class = "btn-success w-100"),
        hr(),
        tags$p(class = "sidebar-title", "Templates"),
        downloadButton("download_template_calib", "Calibration template",
                       class = "btn-sm btn-outline-primary w-100 mb-2"),
        downloadButton("download_template_sample", "Sample template",
                       class = "btn-sm btn-outline-primary w-100")
      ),
      uiOutput("range_validation_alert"),
      layout_column_wrap(
        width = "180px", fill = FALSE,
        value_box("Samples", textOutput("vb_n_samples"), theme = "primary"),
        value_box("In range", textOutput("vb_in_range"), theme = "success"),
        value_box("Below range", textOutput("vb_below"), theme = "warning"),
        value_box("Above range", textOutput("vb_above"), theme = "danger")
      ),
      card(full_screen = TRUE,
           card_header("Calculated concentrations"),
           DTOutput("results_table")),
      card(card_header("Sample data (editable: double-click a cell)"),
           DTOutput("sample_table"))
    )
  ),

  # ---------------------------------------------------------------- Tab 3
  nav_panel(
    "3. Calibration diagnostics",
    two_cols(
      card(card_header("Calibration quality metrics"), uiOutput("calib_qc_metrics")),
      plot_card("Standardized residuals", "residual_plot",
                footer = span(class = "note",
                              "Residuals should scatter randomly within ±2. A curved pattern suggests a",
                              "non-linear response; a funnel shape suggests weighting is needed."))
    ),
    two_cols(
      plot_card("Cook's distance (influence)", "cooks_distance_plot",
                footer = span(class = "note",
                              "Points above 4/n (red line) strongly influence the fit. They are flagged,",
                              "not removed: check them before excluding anything.")),
      plot_card("Back-calculated accuracy", "accuracy_plot",
                footer = span(class = "note",
                              "Acceptance: 85-115% (dashed), 80-120% at the lowest level (dotted)."))
    ),
    card(full_screen = TRUE,
         card_header("Calibration points with diagnostics"),
         DTOutput("calib_diagnostics_table"))
  ),

  # ---------------------------------------------------------------- Tab 4
  nav_panel(
    "4. Sample diagnostics",
    two_cols(
      card(card_header("Sample quality metrics"), uiOutput("sample_qc_metrics")),
      plot_card("Sample responses vs calibration range", "range_coverage_plot",
                footer = span(class = "note",
                              "Grey ticks = calibration standards used; green band = working range.",
                              "Samples outside should be diluted/concentrated or the curve extended."))
    ),
    two_cols(
      plot_card("Concentration distribution", "concentration_distribution"),
      card(full_screen = TRUE, min_height = "360px",
           card_header("Sample accuracy (needs Theoretical_Concentration)"),
           uiOutput("sample_accuracy_info"),
           plotOutput("sample_recovery_plot", height = "320px"))
    ),
    card(card_header("Flagged samples"),
         DTOutput("sample_outliers_table"),
         span(class = "note",
              "Criteria: outside the working range, robust z-score > 3.5 within the group",
              "(median/MAD, see Statistics tab for the grouping column), recovery outside 85-115%."))
  ),

  # ---------------------------------------------------------------- Tab 5
  nav_panel(
    "5. Statistics",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340, title = "Plot and test options",
        selectInput("category", "Grouping column", choices = "Class"),
        selectInput("plot_type", "Plot type",
                    choices = c("Boxplot" = "box", "Violin" = "violin",
                                "Bar (mean ± SE)" = "bar")),
        checkboxInput("show_points", "Show individual points", TRUE),
        checkboxInput("log_scale", "Log scale (y axis)", FALSE),
        checkboxInput("show_stats", "Show mean and n", TRUE),
        checkboxInput("stats_in_range_only", "Use only samples within the working range", FALSE),
        hr(),
        selectizeInput("samples_to_exclude", "Exclude samples", choices = NULL,
                       multiple = TRUE,
                       options = list(placeholder = "Click to choose samples...",
                                      plugins = list("remove_button"))),
        actionButton("reset_exclusions", "Reset exclusions",
                     class = "btn-sm btn-outline-warning w-100"),
        hr(),
        downloadButton("download_plot", "Download plot (PDF)", class = "btn-info w-100")
      ),
      plot_card("Concentration by group", "boxplot", height = "58vh", min_height = "420px"),
      two_cols(
        card(card_header("Group summary"), tableOutput("group_summary")),
        card(card_header("Statistical test"),
             uiOutput("statistical_test"),
             tableOutput("posthoc_table"))
      )
    )
  ),

  # ---------------------------------------------------------------- Guide
  nav_panel(
    "Guide",
    card(
      card_header("Workflow"),
      tags$ol(
        tags$li(tags$b("Calibration:"), " upload the standards, choose model and weighting, check the curve."),
        tags$li(tags$b("Samples:"), " upload sample responses; concentrations are computed immediately and",
                " update automatically whenever the calibration changes."),
        tags$li(tags$b("Diagnostics:"), " check residuals, influential points and back-calculated accuracy",
                " (tab 3), then the sample range coverage and flags (tab 4)."),
        tags$li(tags$b("Statistics:"), " compare groups; exclude samples if justified (and report it).")
      ),
      p(class = "note", "Formula (linear): Concentration = (Response − intercept) / slope × dilution factor.")
    ),
    two_cols(
      card(
        card_header("Weighting and point selection"),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(tags$th("Weighting"), tags$th("When to use"))),
          tags$tbody(
            tags$tr(tags$td(tags$b("None (OLS)")), tags$td("Constant variance (narrow range)")),
            tags$tr(tags$td(tags$b("1/X")), tags$td("Variance increases with concentration")),
            tags$tr(tags$td(tags$b("1/X²")), tags$td("Constant relative error (wide range, LC-MS); usually best over several decades"))
          )
        ),
        tags$ul(
          tags$li(tags$b("Use all points:"), " every standard in the regression."),
          tags$li(tags$b("Bracket the samples:"), " keeps the standards covering the sample responses plus one",
                  " standard below and one above, so no sample is extrapolated."),
          tags$li(tags$b("Manual:"), " tick the points to use."),
          tags$li("Influential points (Cook's distance > 4/n) are only ", tags$b("flagged"), ".",
                  " Removing standards should be a justified decision, not an automatic one.")
        )
      ),
      card(
        card_header("Sample status"),
        tags$table(class = "table table-sm",
          tags$tr(tags$td(tags$b("✓ In range")), tags$td("Response between the lowest and highest standard used")),
          tags$tr(tags$td(tags$b("⚠ Below range")), tags$td("Below the lowest standard (< LLOQ): concentrate or report '< LLOQ'")),
          tags$tr(tags$td(tags$b("⚠ Above range")), tags$td("Above the highest standard: dilute and re-inject")),
          tags$tr(tags$td(tags$b("✗ Not calculable")), tags$td("Response cannot be inverted by the model (e.g. beyond the quadratic maximum)"))
        ),
        p(class = "note", "The 95% CI column (OLS linear model) gives the uncertainty from the calibration",
          " for a single injection.")
      )
    ),
    two_cols(
      card(
        card_header("Calibration metrics"),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(tags$th("Metric"), tags$th("Target"), tags$th("Meaning"))),
          tags$tbody(
            tags$tr(tags$td("R²"), tags$td("≥ 0.99"), tags$td("Necessary but not sufficient: always check residuals and accuracy")),
            tags$tr(tags$td("Back-calculated accuracy"), tags$td("85-115% (80-120% LLOQ)"), tags$td("The most informative check of a calibration")),
            tags$tr(tags$td("Mandel test p"), tags$td("> 0.05"), tags$td("p < 0.05: the quadratic term is significant (curvature)")),
            tags$tr(tags$td("LOD / LOQ"), tags$td("-"), tags$td("3.3·s/slope and 10·s/slope, s = residual SD at the lowest standard"))
          )
        ),
        p(class = "note", "For weighted fits R² is the weighted R²; for the power model it is computed on log-log values.",
          " These are not directly comparable with an unweighted R².")
      ),
      card(
        card_header("Statistics"),
        tags$table(class = "table table-sm table-bordered",
          tags$thead(tags$tr(tags$th("Groups"), tags$th("Parametric"), tags$th("Non-parametric"))),
          tags$tbody(
            tags$tr(tags$td("2"), tags$td("Welch t-test"), tags$td("Mann-Whitney (Wilcoxon)")),
            tags$tr(tags$td("> 2"), tags$td("One-way ANOVA + Tukey HSD"), tags$td("Kruskal-Wallis"))
          )
        ),
        p(class = "note", "With small groups (n < 10) or skewed data, prefer the non-parametric test",
          " or the log scale. *** p < 0.001, ** p < 0.01, * p < 0.05, ns: not significant.")
      )
    ),
    card(
      card_header("Common issues"),
      tags$table(class = "table table-bordered",
        tags$thead(tags$tr(tags$th("Problem"), tags$th("Possible cause"), tags$th("Solution"))),
        tags$tbody(
          tags$tr(tags$td("Response decreases at high levels"), tags$td("Detector saturation"), tags$td("Exclude the saturated levels (manual selection) or dilute")),
          tags$tr(tags$td("Poor accuracy at low levels"), tags$td("High standards dominate the OLS fit"), tags$td("Use 1/X or 1/X² weighting")),
          tags$tr(tags$td("Curved residuals / Mandel p < 0.05"), tags$td("Non-linear response"), tags$td("Narrow the range or use a quadratic model")),
          tags$tr(tags$td("Many samples out of range"), tags$td("Calibration range not adapted"), tags$td("Dilute samples or extend the calibration")),
          tags$tr(tags$td("Large intercept"), tags$td("Blank contamination / background"), tags$td("Check the blank, subtract background"))
        )
      ),
      p(class = "note", "Acceptance criteria are general guidelines (FDA/EMA bioanalytical, ICH Q2):",
        " always refer to the requirements of your field.")
    )
  ),
  nav_spacer(),
  nav_item(tags$span("v3.5", style = "opacity: 0.7;"))
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {

  unit_lab <- reactive({
    u <- trimws(input$conc_unit)
    if (nzchar(u)) paste0(" (", u, ")") else ""
  })
  unit_txt <- reactive(trimws(input$conc_unit))

  # Value boxes show a dash (not the long validation message) until data is loaded
  vb_render <- function(expr, empty = "–", env = parent.frame()) {
    e <- substitute(expr)
    force(env)
    renderText(tryCatch(eval(e, env), error = function(err) empty))
  }

  # =========================================================== CALIBRATION DATA
  calib_raw <- reactiveVal(data.frame(Concentration = rep(NA_real_, 5),
                                      Response = rep(NA_real_, 5)))

  calib_clean <- reactive({
    d <- calib_raw()
    d <- d[is.finite(d$Concentration) & is.finite(d$Response), c("Concentration", "Response")]
    d <- d[order(d$Concentration), , drop = FALSE]
    rownames(d) <- NULL
    d
  })

  load_calib_file <- function(path) {
    df <- tryCatch(detect_and_read_csv(path), error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) {
      showNotification("Could not read the calibration file (need at least 2 columns).", type = "error")
      return()
    }
    nm <- names(df)
    conc_col <- grep("concentration|conc|dose|actual|amount|level", nm, ignore.case = TRUE)[1]
    resp_col <- grep("response|area|signal|intensity|height|absorb", nm, ignore.case = TRUE)[1]
    num_cols <- which(vapply(df, function(x) any(is.finite(to_numeric(x))), logical(1)))
    if (is.na(conc_col)) conc_col <- setdiff(num_cols, resp_col)[1]
    if (is.na(resp_col)) resp_col <- setdiff(num_cols, conc_col)[1]
    if (is.na(conc_col) || is.na(resp_col) || conc_col == resp_col) {
      showNotification("Could not identify the Concentration and Response columns.", type = "error")
      return()
    }
    d <- data.frame(Concentration = to_numeric(df[[conc_col]]),
                    Response = to_numeric(df[[resp_col]]))
    n_bad <- sum(!complete.cases(d))
    calib_raw(d[complete.cases(d), ])

    unit <- regmatches(nm[conc_col], regexpr("\\(([^)]*)\\)", nm[conc_col]))
    if (length(unit)) updateTextInput(session, "conc_unit", value = gsub("[()]", "", unit))

    showNotification(
      sprintf("Loaded %d calibration points: concentration = '%s', response = '%s' (separator: %s)%s",
              sum(complete.cases(d)), nm[conc_col], nm[resp_col], attr(df, "separator"),
              if (n_bad) sprintf(". %d incomplete row(s) ignored.", n_bad) else ""),
      type = "message", duration = 8)
  }
  observeEvent(input$calib_file, load_calib_file(input$calib_file$datapath))

  output$calib_table <- renderDT({
    datatable(calib_raw(), editable = "cell", rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE,
                             scrollY = "300px", scrollCollapse = TRUE))
  })

  observeEvent(input$calib_table_cell_edit, {
    calib_raw(apply_edits(calib_raw(), input$calib_table_cell_edit,
                          c("Concentration", "Response")))
  })
  observeEvent(input$add_calib_row, calib_raw(add_empty_row(calib_raw())))
  observeEvent(input$remove_calib_row, {
    d <- calib_raw()
    if (nrow(d) > 1) calib_raw(d[-nrow(d), , drop = FALSE])
  })

  output$manual_point_selector <- renderUI({
    df <- calib_clean()
    if (nrow(df) == 0) return(helpText("No calibration data yet."))
    checkboxGroupInput(
      "manual_selected_points", "Points to include:",
      choices = setNames(seq_len(nrow(df)),
                         sprintf("C = %s  |  y = %s", fmt(df$Concentration), fmt(df$Response))),
      selected = seq_len(nrow(df))
    )
  })

  # =========================================================== SAMPLE DATA
  samples_raw <- reactiveVal(data.frame(Sample_Name = rep("", 5),
                                        Response = rep(NA_real_, 5),
                                        Class = rep("", 5),
                                        stringsAsFactors = FALSE))
  categories <- reactiveVal("Class")

  samples_clean <- reactive({
    d <- samples_raw()
    d <- d[is.finite(d$Response), , drop = FALSE]
    if (!nrow(d)) return(d)
    d$Sample_Name <- trimws(d$Sample_Name)
    empty <- !nzchar(d$Sample_Name) | is.na(d$Sample_Name)
    d$Sample_Name[empty] <- paste0("Sample_", which(empty))
    d$Sample_Name <- make.unique(d$Sample_Name, sep = "_")
    for (cl in intersect(categories(), names(d))) {
      v <- d[[cl]]
      v[is.na(v) | !nzchar(v)] <- "(none)"
      d[[cl]] <- v
    }
    rownames(d) <- NULL
    d
  })

  load_sample_file <- function(path) {
    df <- tryCatch(detect_and_read_csv(path), error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) {
      showNotification("Could not read the sample file (need at least 2 columns).", type = "error")
      return()
    }
    nm <- names(df)
    name_col <- grep("sample.?name|^name|^sample|^id$|echan|identifier", nm, ignore.case = TRUE)[1]
    if (is.na(name_col)) name_col <- 1
    resp_col <- grep("response|area|signal|intensity|height|absorb", nm, ignore.case = TRUE)[1]
    if (is.na(resp_col)) {
      num_cols <- which(vapply(df, function(x) any(is.finite(to_numeric(x))), logical(1)))
      resp_col <- setdiff(num_cols, name_col)[1]
    }
    if (is.na(resp_col)) {
      showNotification("Could not identify the Response column.", type = "error")
      return()
    }
    theo_col <- grep("theoretical|expected|nominal|spiked", nm, ignore.case = TRUE)[1]
    dil_col  <- grep("dilution", nm, ignore.case = TRUE)[1]

    out <- data.frame(Sample_Name = trimws(as.character(df[[name_col]])),
                      Response = to_numeric(df[[resp_col]]),
                      stringsAsFactors = FALSE)
    if (!is.na(theo_col)) out$Theoretical_Concentration <- to_numeric(df[[theo_col]])
    if (!is.na(dil_col))  out$Dilution_Factor <- to_numeric(df[[dil_col]])

    used <- c(name_col, resp_col, theo_col, dil_col)
    cat_idx <- setdiff(seq_along(nm), used[!is.na(used)])
    cat_idx <- cat_idx[vapply(df[cat_idx], function(x) {
      is.character(x) || is.factor(x) || length(unique(na.omit(x))) <= 20
    }, logical(1))]
    reserved <- c(names(out), "Concentration", "Concentration_Measured", "CI95",
                  "Status", "In_Range", "Recovery_Percent", "Group")
    cats <- character(0)
    for (j in cat_idx) {
      cname <- if (nm[j] %in% reserved) paste0("Cat_", nm[j]) else nm[j]
      out[[cname]] <- as.character(df[[j]])
      cats <- c(cats, cname)
    }
    if (!length(cats)) {
      out$Class <- "All"
      cats <- "Class"
    }

    dups <- sum(duplicated(out$Sample_Name))
    samples_raw(out)
    categories(cats)
    showNotification(
      sprintf("Loaded %d samples (separator: %s). Grouping column(s): %s%s%s%s",
              sum(is.finite(out$Response)), attr(df, "separator"), paste(cats, collapse = ", "),
              if (!is.na(theo_col)) "; theoretical concentrations found" else "",
              if (!is.na(dil_col)) "; dilution factors found" else "",
              if (dups) sprintf(". Warning: %d duplicated name(s) were made unique.", dups) else ""),
      type = if (dups) "warning" else "message", duration = 8)
  }
  observeEvent(input$sample_file, load_sample_file(input$sample_file$datapath))

  # Example data shipped next to the app (CalibUV2.csv / SampleUV2.csv)
  example_files <- c(calib = "CalibUV2.csv", sample = "SampleUV2.csv")
  output$example_ui <- renderUI({
    if (!all(file.exists(example_files))) return(NULL)
    actionButton("load_example", "Load example data", class = "btn-sm btn-outline-secondary w-100")
  })
  observeEvent(input$load_example, {
    load_calib_file(example_files[["calib"]])
    load_sample_file(example_files[["sample"]])
  })

  output$sample_table <- renderDT({
    datatable(samples_raw(), editable = "cell", rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE,
                             scrollY = "300px", scrollX = TRUE, scrollCollapse = TRUE))
  })
  observeEvent(input$sample_table_cell_edit, {
    samples_raw(apply_edits(samples_raw(), input$sample_table_cell_edit, SAMPLE_NUMERIC))
  })
  observeEvent(input$add_sample_row, samples_raw(add_empty_row(samples_raw())))
  observeEvent(input$remove_sample_row, {
    d <- samples_raw()
    if (nrow(d) > 1) samples_raw(d[-nrow(d), , drop = FALSE])
  })

  # =========================================================== FIT (reactive)
  selection <- reactive({
    df <- calib_clean()
    validate(need(nrow(df) >= 3, "Load or enter at least 3 calibration points (tab 1)."))
    mode <- input$point_selection_mode
    if (mode == "all") {
      list(used = rep(TRUE, nrow(df)), note = "All calibration points are used.")
    } else if (mode == "manual") {
      idx <- as.integer(input$manual_selected_points)
      used <- seq_len(nrow(df)) %in% idx
      validate(need(sum(used) >= 3, "Manual selection: tick at least 3 points."))
      list(used = used, note = sprintf("%d point(s) selected manually.", sum(used)))
    } else {
      select_intelligent(df, samples_clean()$Response, input$safety_margin / 100,
                         input$include_zero, input$min_points)
    }
  })

  fit <- reactive({
    df <- calib_clean()
    sel <- selection()
    d <- df[sel$used, , drop = FALSE]

    f <- tryCatch(fit_calibration(d, input$regression_type, input$weighting_scheme),
                  error = function(e) e)
    validate(need(!inherits(f, "error"), if (inherits(f, "error")) conditionMessage(f)))

    cd <- cooks.distance(f$model)
    thr <- 4 / length(cd)
    flagged <- which(cd > thr)
    removed <- NULL
    if (isTRUE(input$drop_influential) && length(flagged)) {
      removed <- f$data[flagged, c("Concentration", "Response")]
      f2 <- tryCatch(fit_calibration(f$data[-flagged, ], input$regression_type,
                                     input$weighting_scheme),
                     error = function(e) e)
      validate(need(!inherits(f2, "error"),
                    "Too few points left after excluding the influential points."))
      f2$notes <- unique(c(f$notes, f2$notes))
      f <- f2
      cd <- cooks.distance(f$model)
      thr <- 4 / length(cd)
    }

    pos <- f$data$Concentration > 0
    f$cooks <- unname(cd)
    f$cooks_thr <- thr
    f$removed <- removed
    f$all <- df
    f$sel_note <- sel$note
    f$range_resp <- range(f$data$Response[pos])
    f$range_conc <- range(f$data$Concentration[pos])
    f$r2 <- summary(f$model)$r.squared
    f$r2_label <- if (f$type == "power") "R² (log-log)" else if (f$weighted) "R² (weighted)" else "R²"
    f$lod <- lod_loq(f)
    f$mandel <- mandel_p(f)
    bc <- conc_from_response(f, f$data$Response)
    f$data$Back_Calc <- bc
    f$data$Accuracy <- ifelse(pos, 100 * bc / f$data$Concentration, NA_real_)
    f
  })

  saturation_warning <- reactive({
    df <- calib_clean()
    if (nrow(df) < 3) return(NULL)
    a <- aggregate(Response ~ Concentration, data = df, FUN = mean)
    dec <- which(diff(a$Response) <= 0)
    if (!length(dec)) return(NULL)
    sprintf(paste("The response does not increase between C = %s and C = %s%s:",
                  "possible detector saturation or preparation error.",
                  "Do not use these levels for quantification (manual selection) unless justified."),
            fmt(a$Concentration[dec[1]]), fmt(a$Concentration[dec[1] + 1]), unit_lab())
  })

  weight_label <- function(f) {
    switch(f$scheme, none = "unweighted (OLS)", "1/x" = "1/X weighted", "1/x2" = "1/X² weighted")
  }

  # =========================================================== TAB 1 OUTPUTS
  output$calib_warnings <- renderUI({
    msgs <- list()
    sw <- saturation_warning()
    if (!is.null(sw)) msgs <- c(msgs, list(div(class = "alert alert-warning py-2 mb-2", tags$b("⚠ "), sw)))
    f <- tryCatch(fit(), error = function(e) NULL)
    if (!is.null(f)) {
      for (n in f$notes) msgs <- c(msgs, list(div(class = "alert alert-info py-2 mb-2", n)))
      bad <- sum(f$data$Accuracy < 85 | f$data$Accuracy > 115, na.rm = TRUE)
      if (bad) {
        msgs <- c(msgs, list(div(class = "alert alert-warning py-2 mb-2",
          sprintf("⚠ %d standard(s) back-calculate outside 85-115%% (see tab 3). Consider weighting, another model or point selection.", bad))))
      }
    }
    tagList(msgs)
  })

  output$vb_r2_label <- vb_render(fit()$r2_label, empty = "R²")
  output$vb_r2 <- vb_render(sprintf("%.4f", fit()$r2))
  output$vb_points <- vb_render(sprintf("%d / %d", nrow(fit()$data), nrow(fit()$all)))
  output$vb_range <- vb_render({
    f <- fit()
    paste(fmt(f$range_conc[1], 3), "–", fmt(f$range_conc[2], 3), unit_txt())
  })
  output$vb_lod <- vb_render({
    l <- fit()$lod
    if (is.null(l)) "n/a" else paste(fmt(l$lod, 3), "/", fmt(l$loq, 3), unit_txt())
  })

  calib_plot_obj <- reactive({
    f <- fit()
    logax <- isTRUE(input$calib_log_axes)
    all <- f$all
    key <- function(d) paste(d$Concentration, d$Response)
    all$Status <- ifelse(key(all) %in% key(f$data), "Used",
                  ifelse(key(all) %in% key(f$removed), "Influential (excluded)", "Not used"))
    pts <- if (isTRUE(input$show_excluded_points)) all else all[all$Status == "Used", ]

    xr <- range(f$data$Concentration)
    pos_min <- min(f$data$Concentration[f$data$Concentration > 0])
    xs <- if (logax || f$type %in% c("log", "power")) {
      10^seq(log10(pos_min), log10(xr[2]), length.out = 300)
    } else {
      seq(xr[1], xr[2], length.out = 300)
    }
    pred <- data.frame(Concentration = xs, Response = predict_response(f, xs))

    smp <- NULL
    if (isTRUE(input$show_samples_on_curve)) {
      r <- tryCatch(results(), error = function(e) NULL)
      if (!is.null(r)) {
        smp <- r[is.finite(r$Concentration_Measured), ]
        smp <- data.frame(Concentration = smp$Concentration_Measured,
                          Response = smp$Response, Status = "Samples")
      }
    }
    if (logax) {
      keep_pos <- function(d) if (is.null(d)) NULL else d[d$Concentration > 0 & d$Response > 0, ]
      pts <- keep_pos(pts); pred <- keep_pos(pred); smp <- keep_pos(smp)
    }

    band_x <- if (logax) range(c(pts$Concentration, pred$Concentration)) else c(-Inf, Inf)
    p <- ggplot() +
      annotate("rect", xmin = band_x[1], xmax = band_x[2],
               ymin = f$range_resp[1], ymax = f$range_resp[2],
               fill = "#27ae60", alpha = 0.06) +
      geom_line(data = pred, aes(Concentration, Response), colour = "#27ae60", linewidth = 1.3) +
      geom_point(data = pts, aes(Concentration, Response, colour = Status, shape = Status),
                 size = 3.5, stroke = 1.5)
    if (!is.null(smp) && nrow(smp)) {
      p <- p + geom_point(data = smp, aes(Concentration, Response, colour = Status, shape = Status),
                          size = 3, alpha = 0.8)
    }
    p <- p +
      scale_colour_manual(values = c("Used" = "#2C3E50", "Not used" = "grey65",
                                     "Influential (excluded)" = "#e74c3c", "Samples" = "#e67e22"),
                          name = NULL) +
      scale_shape_manual(values = c("Used" = 16, "Not used" = 1,
                                    "Influential (excluded)" = 4, "Samples" = 17),
                         name = NULL)
    if (logax) p <- p + scale_x_log10() + scale_y_log10()

    p + labs(
      title = sprintf("Calibration curve: %s model, %s", f$type, weight_label(f)),
      subtitle = sprintf("%s = %.4f  |  %d of %d points used  |  green band = working range",
                         f$r2_label, f$r2, nrow(f$data), nrow(f$all)),
      x = paste0("Concentration", unit_lab()), y = "Response",
      caption = if (logax) {
        "Log-log axes: zero values are not shown."
      } else if (diff(log10(f$range_conc)) > 2) {
        "Tip: the range spans more than 2 decades - tick 'Log-log axes' to see the low standards."
      } else NULL
    ) + theme_app()
  })

  output$calib_plot <- renderPlot(calib_plot_obj(), res = 96)

  output$regression_equation <- renderUI({
    f <- fit()
    b <- unname(coef(f$model))
    eq <- switch(f$type,
      linear    = sprintf("y = %s·x %s", fmt(b[2]), signed(b[1])),
      log       = sprintf("y = %s·ln(x) %s", fmt(b[2]), signed(b[1])),
      quadratic = sprintf("y = %s·x² %s·x %s", fmt(b[3]), signed(b[2]), signed(b[1])),
      power     = sprintf("y = %s·x^%s", fmt(exp(b[1])), fmt(b[2]))
    )
    quality <- if (f$r2 >= 0.99) c("Excellent", "#27ae60") else if (f$r2 >= 0.95) c("Acceptable", "#f39c12") else c("Poor", "#e74c3c")
    tagList(
      div(class = "eq", eq),
      p(class = "note", sprintf("y = response, x = concentration%s; %s", unit_lab(), weight_label(f))),
      p(tags$b(sprintf("%s = %.4f", f$r2_label, f$r2)),
        tags$span(sprintf(" (%s)", quality[1]), style = sprintf("color:%s; font-weight:600;", quality[2]))),
      if (is.finite(f$mandel)) {
        p(sprintf("Mandel test (linear vs quadratic): p = %s", fmt_p(f$mandel)),
          if (f$mandel < 0.05) tags$span(" → significant curvature", style = "color:#e67e22; font-weight:600;"))
      }
    )
  })

  output$point_selection_summary <- renderUI({
    f <- fit()
    flagged <- sum(f$cooks > f$cooks_thr, na.rm = TRUE)
    tagList(
      hr(),
      p(f$sel_note),
      tags$table(class = "table table-sm metric-table",
        tags$tr(tags$td("Calibration points available"), tags$td(tags$b(nrow(f$all)))),
        tags$tr(tags$td("Points used in the fit"), tags$td(tags$b(nrow(f$data)))),
        tags$tr(tags$td("Influential points flagged (Cook's D > 4/n)"),
                tags$td(tags$b(flagged, style = if (flagged) "color:#e74c3c;" else ""))),
        if (!is.null(f$removed)) {
          tags$tr(tags$td("Points automatically excluded"),
                  tags$td(tags$b(nrow(f$removed), style = "color:#e74c3c;"),
                          sprintf(" (C = %s)", paste(fmt(f$removed$Concentration), collapse = ", "))))
        }
      )
    )
  })

  # =========================================================== RESULTS (reactive)
  results <- reactive({
    f <- fit()
    s <- samples_clean()
    validate(need(nrow(s) > 0, "Load or enter sample data (tab 2)."))

    conc <- conc_from_response(f, s$Response)
    dil <- if ("Dilution_Factor" %in% names(s)) {
      ifelse(is.finite(s$Dilution_Factor) & s$Dilution_Factor > 0, s$Dilution_Factor, 1)
    } else 1
    s$Concentration_Measured <- conc
    s$Concentration <- conc * dil
    half <- ci_halfwidth(f, s$Response)
    if (!is.null(half)) s$CI95 <- half * dil

    s$Status <- ifelse(!is.finite(conc), STATUS_LEVELS[4],
                ifelse(s$Response < f$range_resp[1], STATUS_LEVELS[2],
                ifelse(s$Response > f$range_resp[2], STATUS_LEVELS[3], STATUS_LEVELS[1])))
    s$In_Range <- s$Status == STATUS_LEVELS[1]

    if ("Theoretical_Concentration" %in% names(s)) {
      th <- s$Theoretical_Concentration
      s$Recovery_Percent <- ifelse(is.finite(th) & th > 0, 100 * s$Concentration / th, NA_real_)
    }
    s
  })

  status_counts <- reactive({
    r <- results()
    setNames(vapply(STATUS_LEVELS, function(l) sum(r$Status == l), numeric(1)), STATUS_LEVELS)
  })

  output$vb_n_samples <- vb_render(nrow(results()))
  output$vb_in_range <- vb_render(status_counts()[[1]])
  output$vb_below <- vb_render(status_counts()[[2]])
  output$vb_above <- vb_render(status_counts()[[3]] + status_counts()[[4]])

  output$range_validation_alert <- renderUI({
    r <- results()
    n_out <- sum(!r$In_Range)
    if (n_out == 0) {
      div(class = "alert alert-success py-2", tags$b("✓ All samples are within the working range."))
    } else {
      div(class = "alert alert-warning py-2",
          tags$b(sprintf("⚠ %d of %d samples are outside the working range.", n_out, nrow(r))),
          " Above range: dilute and re-inject. Below range: report as < LLOQ or concentrate.",
          " In 'bracket the samples' mode the curve adapts automatically to the sample responses.")
    }
  })

  display_results <- reactive({
    r <- results()
    cats <- intersect(categories(), names(r))
    cols <- c("Sample_Name", cats, "Response",
              if ("Dilution_Factor" %in% names(r)) c("Dilution_Factor", "Concentration_Measured"),
              "Concentration", if ("CI95" %in% names(r)) "CI95", "Status",
              if ("Theoretical_Concentration" %in% names(r)) c("Theoretical_Concentration", "Recovery_Percent"))
    d <- r[, cols, drop = FALSE]
    u <- unit_lab()
    rename <- c(Sample_Name = "Sample", Dilution_Factor = "Dilution factor",
                Concentration_Measured = paste0("Conc. in vial", u),
                Concentration = paste0("Concentration", u),
                CI95 = "95% CI (±)", Theoretical_Concentration = paste0("Theoretical", u),
                Recovery_Percent = "Recovery (%)")
    hit <- names(d) %in% names(rename)
    names(d)[hit] <- rename[names(d)[hit]]
    d
  })

  output$results_table <- renderDT({
    d <- display_results()
    num <- which(vapply(d, is.numeric, logical(1)))
    dt <- datatable(d, rownames = FALSE, selection = "none",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
                                   buttons = c("copy"))) %>%
      formatSignif(columns = num, digits = 4) %>%
      formatStyle("Status", backgroundColor = styleEqual(STATUS_LEVELS, STATUS_COLOURS))
    if ("Recovery (%)" %in% names(d)) {
      dt <- dt %>% formatRound("Recovery (%)", 1) %>%
        formatStyle("Recovery (%)", backgroundColor = styleInterval(
          c(85, 115), c("#f8d7da", "#d4edda", "#f8d7da")))
    }
    dt
  })

  # =========================================================== TAB 3 OUTPUTS
  output$calib_qc_metrics <- renderUI({
    f <- fit()
    acc <- f$data$Accuracy
    pos <- f$data$Concentration > 0
    lowest <- which(pos)[which.min(f$data$Concentration[pos])]
    n_bad <- sum((acc < 85 | acc > 115)[-lowest], na.rm = TRUE) +
      sum(acc[lowest] < 80 | acc[lowest] > 120, na.rm = TRUE)
    row <- function(label, value, style = "") tags$tr(tags$td(label), tags$td(tags$b(value, style = style)))
    r2_col <- if (f$r2 >= 0.99) "color:#27ae60;" else if (f$r2 >= 0.95) "color:#f39c12;" else "color:#e74c3c;"
    tagList(
      tags$table(class = "table table-sm metric-table",
        row("Model", tools::toTitleCase(f$type)),
        row("Weighting", weight_label(f)),
        row(f$r2_label, sprintf("%.4f", f$r2), r2_col),
        row("Points used", nrow(f$data)),
        row("Working range (concentration)",
            paste(fmt(f$range_conc[1]), "–", fmt(f$range_conc[2]), unit_txt())),
        row("Working range (response)", paste(fmt(f$range_resp[1]), "–", fmt(f$range_resp[2]))),
        if (!f$weighted && f$type != "power") row("Residual SD, s(y/x)", fmt(summary(f$model)$sigma)),
        row("Mean |bias| of back-calculated standards", sprintf("%.1f %%", mean(abs(acc - 100), na.rm = TRUE))),
        row("Standards outside acceptance", sprintf("%d / %d", n_bad, sum(pos)),
            if (n_bad) "color:#e74c3c;" else "color:#27ae60;"),
        row("Mandel test (quadratic term) p", if (is.finite(f$mandel)) fmt_p(f$mandel) else "n/a (linear model only)"),
        row("LOD / LOQ", if (is.null(f$lod)) "n/a (linear/quadratic only)" else
              paste(fmt(f$lod$lod, 3), "/", fmt(f$lod$loq, 3), unit_txt()))
      ),
      if (n_bad == 0 && f$r2 >= 0.99) {
        div(class = "alert alert-success py-2", "✓ Calibration meets the usual acceptance criteria.")
      } else {
        div(class = "alert alert-warning py-2",
            "⚠ Check the flagged standards: try 1/X² weighting, another model, or a narrower range.")
      }
    )
  })

  output$residual_plot <- renderPlot({
    f <- fit()
    d <- f$data
    d$Std_Res <- rstandard(f$model)
    ggplot(d, aes(Concentration, Std_Res)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2, fill = "#27ae60", alpha = 0.07) +
      geom_hline(yintercept = 0, colour = "#e74c3c", linetype = "dashed", linewidth = 0.9) +
      geom_hline(yintercept = c(-2, 2), colour = "grey60", linetype = "dotted") +
      geom_point(size = 3.5, colour = "#2C3E50") +
      conc_x_scale(d$Concentration, name = paste0("Concentration", unit_lab())) +
      labs(title = "Standardized residuals", y = "Standardized residual",
           subtitle = weight_label(f)) +
      theme_app()
  }, res = 96)

  output$cooks_distance_plot <- renderPlot({
    f <- fit()
    d <- f$data
    d$Cooks <- f$cooks
    d$Label <- factor(fmt(d$Concentration, 3), levels = unique(fmt(d$Concentration, 3)))
    d$Flag <- ifelse(d$Cooks > f$cooks_thr, "Influential", "OK")
    ggplot(d, aes(Label, Cooks, fill = Flag)) +
      geom_col(alpha = 0.85) +
      geom_hline(yintercept = f$cooks_thr, colour = "#e74c3c", linetype = "dashed", linewidth = 1) +
      scale_fill_manual(values = c("Influential" = "#e74c3c", "OK" = "#3498db"), name = NULL) +
      labs(title = "Cook's distance", subtitle = sprintf("Threshold 4/n = %.3f", f$cooks_thr),
           x = paste0("Concentration", unit_lab()), y = "Cook's distance") +
      theme_app() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, res = 96)

  output$accuracy_plot <- renderPlot({
    f <- fit()
    d <- f$data[f$data$Concentration > 0 & is.finite(f$data$Accuracy), ]
    validate(need(nrow(d) > 0, "No non-zero standards to evaluate."))
    ggplot(d, aes(Concentration, Accuracy)) +
      geom_hline(yintercept = 100, colour = "#27ae60", linewidth = 1) +
      geom_hline(yintercept = c(85, 115), colour = "#f39c12", linetype = "dashed") +
      geom_hline(yintercept = c(80, 120), colour = "#f39c12", linetype = "dotted") +
      geom_line(colour = "#3498db", alpha = 0.5) +
      geom_point(aes(colour = Accuracy >= 85 & Accuracy <= 115), size = 3.5) +
      scale_colour_manual(values = c("TRUE" = "#2C3E50", "FALSE" = "#e74c3c"), guide = "none") +
      conc_x_scale(d$Concentration, name = paste0("Nominal concentration", unit_lab())) +
      labs(title = "Back-calculated standards", y = "Accuracy (%)") +
      theme_app()
  }, res = 96)

  output$calib_diagnostics_table <- renderDT({
    f <- fit()
    d <- f$data
    pred <- predict_response(f, d$Concentration)
    out <- data.frame(
      Concentration = d$Concentration,
      Response = d$Response,
      Predicted = pred,
      Residual = d$Response - pred,
      `Residual (%)` = ifelse(d$Response != 0, 100 * (d$Response - pred) / d$Response, NA),
      `Std. residual` = unname(rstandard(f$model)),
      `Back-calculated` = d$Back_Calc,
      `Accuracy (%)` = d$Accuracy,
      Weight = d$.w,
      `Cook's D` = f$cooks,
      check.names = FALSE
    )
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, scrollX = TRUE)) %>%
      formatSignif(columns = c(1:4, 7, 9), digits = 4) %>%
      formatRound(columns = c(5, 6, 8), digits = 1) %>%
      formatRound(columns = 10, digits = 3) %>%
      formatStyle("Cook's D", backgroundColor = styleInterval(f$cooks_thr, c("white", "#f8d7da"))) %>%
      formatStyle("Accuracy (%)", backgroundColor = styleInterval(
        c(80, 85, 115, 120), c("#f8d7da", "#fff3cd", "#d4edda", "#fff3cd", "#f8d7da")))
  })

  # =========================================================== GROUPING
  observe({
    cats <- categories()
    cur <- isolate(input$category)
    updateSelectInput(session, "category", choices = cats,
                      selected = if (!is.null(cur) && cur %in% cats) cur else cats[1])
  })

  group_col <- reactive({
    cats <- categories()
    if (!is.null(input$category) && input$category %in% cats) input$category else cats[1]
  })

  observe({
    r <- tryCatch(results(), error = function(e) NULL)
    nm <- if (is.null(r)) character(0) else r$Sample_Name
    updateSelectizeInput(session, "samples_to_exclude", choices = nm,
                         selected = intersect(isolate(input$samples_to_exclude), nm))
  })
  observeEvent(input$reset_exclusions, {
    updateSelectizeInput(session, "samples_to_exclude", selected = character(0))
  })

  # =========================================================== TAB 4 OUTPUTS
  output$sample_qc_metrics <- renderUI({
    r <- results()
    sc <- status_counts()
    n <- nrow(r)
    pct <- function(k) sprintf("%d / %d (%.0f%%)", as.integer(k), n, 100 * k / n)
    row <- function(label, value, style = "") tags$tr(tags$td(label), tags$td(tags$b(value, style = style)))
    tagList(
      tags$table(class = "table table-sm metric-table",
        row("Number of samples", n),
        row("In range", pct(sc[[1]]), if (sc[[1]] == n) "color:#27ae60;" else "color:#f39c12;"),
        row("Below range (< lowest standard)", pct(sc[[2]])),
        row("Above range (> highest standard)", pct(sc[[3]])),
        row("Not calculable", pct(sc[[4]])),
        if ("Dilution_Factor" %in% names(r)) row("Dilution factors applied", "yes")
      ),
      if ("Recovery_Percent" %in% names(r)) {
        rec <- r$Recovery_Percent
        tagList(
          h6("Accuracy (samples with theoretical concentration)"),
          tags$table(class = "table table-sm metric-table",
            row("Mean recovery", sprintf("%.1f %%", mean(rec, na.rm = TRUE))),
            row("Recovery SD", sprintf("%.1f %%", sd(rec, na.rm = TRUE))),
            row("Within 85-115%", sprintf("%d / %d", sum(rec >= 85 & rec <= 115, na.rm = TRUE), sum(is.finite(rec))))
          )
        )
      }
    )
  })

  output$range_coverage_plot <- renderPlot({
    r <- results()
    f <- fit()
    r <- r[r$Response > 0, ]
    r$Status <- factor(r$Status, levels = STATUS_LEVELS)
    std <- f$data$Response[f$data$Response > 0]
    ggplot(r, aes(Response, 0, colour = Status)) +
      annotate("rect", xmin = f$range_resp[1], xmax = f$range_resp[2], ymin = -Inf, ymax = Inf,
               fill = "#27ae60", alpha = 0.10) +
      geom_vline(xintercept = std, colour = "grey60", linewidth = 0.6) +
      geom_jitter(width = 0, height = 0.3, size = 3.5, alpha = 0.85) +
      scale_x_log10(labels = function(b) fmt(b, 3)) +
      scale_y_continuous(limits = c(-0.6, 0.6), breaks = NULL) +
      scale_colour_manual(values = setNames(c("#27ae60", "#f39c12", "#e74c3c", "grey40"), STATUS_LEVELS),
                          drop = FALSE, name = NULL) +
      labs(title = "Sample responses vs working range", x = "Response (log scale)", y = NULL) +
      theme_app()
  }, res = 96)

  output$concentration_distribution <- renderPlot({
    r <- results()
    g <- group_col()
    r <- r[is.finite(r$Concentration), ]
    validate(need(nrow(r) > 0, "No calculable concentrations."))
    r$Group <- if (g %in% names(r)) r[[g]] else "All"
    ggplot(r, aes(Concentration, fill = Group)) +
      geom_histogram(bins = 15, alpha = 0.65, position = "identity", colour = "white") +
      conc_x_scale(r$Concentration, name = paste0("Concentration", unit_lab())) +
      scale_fill_brewer(palette = "Set2", name = g) +
      labs(title = "Concentration distribution", y = "Count") +
      theme_app()
  }, res = 96)

  output$sample_accuracy_info <- renderUI({
    r <- results()
    if (!"Recovery_Percent" %in% names(r)) {
      return(div(class = "alert alert-secondary py-2",
                 "No theoretical concentrations in the sample file. Add a 'Theoretical_Concentration'",
                 " column (e.g. QC or spiked samples) to assess accuracy."))
    }
    rec <- r$Recovery_Percent
    p(sprintf("Mean recovery %.1f %% (SD %.1f %%), %d / %d within 85-115%%.",
              mean(rec, na.rm = TRUE), sd(rec, na.rm = TRUE),
              sum(rec >= 85 & rec <= 115, na.rm = TRUE), sum(is.finite(rec))))
  })

  output$sample_recovery_plot <- renderPlot({
    r <- results()
    validate(need("Recovery_Percent" %in% names(r), "No theoretical concentrations provided."))
    r <- r[is.finite(r$Recovery_Percent), ]
    validate(need(nrow(r) > 0, "No sample with a theoretical concentration."))
    r$Sample_Name <- factor(r$Sample_Name, levels = r$Sample_Name)
    r$Within <- ifelse(r$Recovery_Percent >= 85 & r$Recovery_Percent <= 115, "Within 85-115%", "Outside limits")
    ggplot(r, aes(Sample_Name, Recovery_Percent, colour = Within)) +
      geom_hline(yintercept = 100, colour = "#27ae60", linewidth = 1) +
      geom_hline(yintercept = c(85, 115), colour = "#f39c12", linetype = "dashed") +
      geom_point(size = 3.5) +
      scale_colour_manual(values = c("Within 85-115%" = "#3498db", "Outside limits" = "#e74c3c"), name = NULL) +
      labs(x = NULL, y = "Recovery (%)") +
      theme_app() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, res = 96)

  output$sample_outliers_table <- renderDT({
    r <- results()
    g <- group_col()
    grp <- if (g %in% names(r)) r[[g]] else rep("All", nrow(r))
    rz <- ave(r$Concentration, grp, FUN = function(x) {
      m <- median(x, na.rm = TRUE)
      s <- mad(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else (x - m) / s
    })
    reasons <- vapply(seq_len(nrow(r)), function(i) {
      why <- character(0)
      if (!r$In_Range[i]) why <- c(why, sub("^\\S+ ", "", r$Status[i]))
      if (is.finite(rz[i]) && abs(rz[i]) > 3.5) why <- c(why, sprintf("robust z = %.1f in %s", rz[i], grp[i]))
      if ("Recovery_Percent" %in% names(r) && is.finite(r$Recovery_Percent[i]) &&
          (r$Recovery_Percent[i] < 85 || r$Recovery_Percent[i] > 115)) {
        why <- c(why, sprintf("recovery %.0f%%", r$Recovery_Percent[i]))
      }
      paste(why, collapse = "; ")
    }, character(1))
    flag <- nzchar(reasons)
    if (!any(flag)) {
      out <- data.frame(Message = "No flagged samples: all samples pass the QC criteria.")
      return(datatable(out, rownames = FALSE, options = list(dom = "t")))
    }
    out <- data.frame(Sample = r$Sample_Name, Group = grp, Response = r$Response,
                      Concentration = r$Concentration, `Robust z` = rz,
                      Reason = reasons, check.names = FALSE)[flag, ]
    datatable(out, rownames = FALSE, selection = "none",
              options = list(dom = "tp", pageLength = 15, scrollX = TRUE)) %>%
      formatSignif(columns = c("Response", "Concentration"), digits = 4) %>%
      formatRound("Robust z", 2)
  })

  # =========================================================== TAB 5 STATISTICS
  stats_data <- reactive({
    r <- results()
    g <- group_col()
    validate(need(g %in% names(r), "Choose a grouping column."))
    d <- r[!r$Sample_Name %in% input$samples_to_exclude, , drop = FALSE]
    if (isTRUE(input$stats_in_range_only)) d <- d[d$In_Range, , drop = FALSE]
    d <- d[is.finite(d$Concentration), , drop = FALSE]
    if (isTRUE(input$log_scale)) d <- d[d$Concentration > 0, , drop = FALSE]
    validate(need(nrow(d) > 0, "No samples left after exclusions/filters."))
    d$Group <- factor(d[[g]])
    d
  })

  stats_plot_obj <- reactive({
    d <- stats_data()
    g <- group_col()
    n_total <- nrow(results())
    if (input$plot_type == "bar") {
      s <- d %>% group_by(Group) %>%
        summarise(Mean = mean(Concentration), SE = sd(Concentration) / sqrt(n()), .groups = "drop") %>%
        mutate(SE = ifelse(is.finite(SE), SE, 0),
               Low = if (isTRUE(input$log_scale)) pmax(Mean - SE, Mean * 0.05) else Mean - SE)
      p <- ggplot(s, aes(Group, Mean, fill = Group)) +
        geom_col(alpha = 0.75, width = 0.6) +
        geom_errorbar(aes(ymin = Low, ymax = Mean + SE), width = 0.2)
    } else {
      p <- ggplot(d, aes(Group, Concentration, fill = Group))
      p <- p + if (input$plot_type == "violin") {
        geom_violin(alpha = 0.65, trim = FALSE)
      } else {
        geom_boxplot(alpha = 0.65, width = 0.6, outlier.shape = if (isTRUE(input$show_points)) NA else 19)
      }
    }
    if (isTRUE(input$show_points)) {
      p <- p + geom_jitter(data = d, aes(Group, Concentration), inherit.aes = FALSE,
                           width = 0.12, height = 0, size = 2.8, alpha = 0.7, colour = "#2C3E50")
    }
    if (isTRUE(input$show_stats)) {
      lab <- d %>% group_by(Group) %>%
        summarise(Mean = mean(Concentration), Top = max(Concentration), n = n(), .groups = "drop")
      p <- p + geom_text(data = lab, aes(Group, Top, label = sprintf("mean %s\nn = %d", fmt(Mean, 3), n)),
                         inherit.aes = FALSE, vjust = -0.4, size = 4.5, fontface = "bold")
    }
    yexp <- expansion(mult = c(0.05, if (isTRUE(input$show_stats)) 0.3 else 0.08))
    p <- p + if (isTRUE(input$log_scale)) scale_y_log10(expand = yexp) else scale_y_continuous(expand = yexp)
    n_exc <- n_total - nrow(d)
    p + scale_fill_brewer(palette = "Set2", guide = "none") +
      labs(title = sprintf("Concentration by %s", g),
           subtitle = sprintf("n = %d samples%s", nrow(d),
                              if (n_exc > 0) sprintf(" (%d excluded or filtered)", n_exc) else ""),
           x = g, y = paste0("Concentration", unit_lab(), if (isTRUE(input$log_scale)) ", log scale" else "")) +
      theme_app(16)
  })

  output$boxplot <- renderPlot(stats_plot_obj(), res = 96)

  output$group_summary <- renderTable({
    d <- stats_data()
    d %>% group_by(Group) %>%
      summarise(n = n(),
                Mean = fmt(mean(Concentration)),
                SD = fmt(sd(Concentration)),
                Median = fmt(median(Concentration)),
                `CV (%)` = sprintf("%.1f", 100 * sd(Concentration) / mean(Concentration)),
                .groups = "drop") %>%
      as.data.frame()
  }, striped = TRUE, hover = TRUE, width = "100%")

  test_info <- reactive({
    d <- droplevels(stats_data())
    k <- nlevels(d$Group)
    ns <- table(d$Group)
    if (k < 2) return(list(msg = "At least 2 groups are needed for a statistical comparison."))
    if (any(ns < 2)) return(list(msg = sprintf("Each group needs at least 2 samples (%s).",
                                               paste(names(ns), ns, sep = ": n = ", collapse = ", "))))
    if (k == 2) {
      tt <- t.test(Concentration ~ Group, data = d)
      wt <- suppressWarnings(wilcox.test(Concentration ~ Group, data = d, exact = FALSE))
      x1 <- d$Concentration[d$Group == levels(d$Group)[1]]
      x2 <- d$Concentration[d$Group == levels(d$Group)[2]]
      sp <- sqrt(((length(x1) - 1) * var(x1) + (length(x2) - 1) * var(x2)) / (length(x1) + length(x2) - 2))
      list(k = 2, rows = list(
        c("Welch t-test", sprintf("t = %.3f, df = %.1f", tt$statistic, tt$parameter), fmt_p(tt$p.value), p_stars(tt$p.value)),
        c("Mann-Whitney (Wilcoxon)", sprintf("W = %.0f", wt$statistic), fmt_p(wt$p.value), p_stars(wt$p.value)),
        c("Cohen's d (pooled SD)", sprintf("%.2f", abs(mean(x1) - mean(x2)) / sp), "", "")
      ))
    } else {
      a <- aov(Concentration ~ Group, data = d)
      sa <- summary(a)[[1]]
      kw <- kruskal.test(Concentration ~ Group, data = d)
      tk <- TukeyHSD(a)$Group
      list(k = k, rows = list(
        c("One-way ANOVA", sprintf("F(%d, %d) = %.3f", sa$Df[1], sa$Df[2], sa$`F value`[1]),
          fmt_p(sa$`Pr(>F)`[1]), p_stars(sa$`Pr(>F)`[1])),
        c("Kruskal-Wallis", sprintf("χ² = %.3f, df = %d", kw$statistic, kw$parameter),
          fmt_p(kw$p.value), p_stars(kw$p.value))
      ), tukey = data.frame(
        Comparison = rownames(tk),
        Difference = fmt(tk[, "diff"]),
        `95% CI` = paste(fmt(tk[, "lwr"]), "to", fmt(tk[, "upr"])),
        `p adj` = vapply(tk[, "p adj"], fmt_p, character(1)),
        ` ` = vapply(tk[, "p adj"], p_stars, character(1)),
        check.names = FALSE
      ))
    }
  })

  output$statistical_test <- renderUI({
    ti <- test_info()
    if (!is.null(ti$msg)) return(div(class = "alert alert-secondary py-2", ti$msg))
    n_exc <- length(input$samples_to_exclude)
    tagList(
      if (n_exc > 0) div(class = "alert alert-warning py-1", sprintf("%d sample(s) excluded manually: report this.", n_exc)),
      tags$table(class = "table table-sm table-striped",
        tags$thead(tags$tr(tags$th("Test"), tags$th("Statistic"), tags$th("p-value"), tags$th(""))),
        tags$tbody(lapply(ti$rows, function(r) tags$tr(lapply(r, tags$td))))
      ),
      if (ti$k > 2) h6("Tukey HSD post-hoc comparisons")
    )
  })

  output$posthoc_table <- renderTable({
    ti <- test_info()
    if (is.null(ti$tukey)) return(NULL)
    ti$tukey
  }, striped = TRUE, width = "100%")

  # =========================================================== DOWNLOADS
  output$download_template_calib <- downloadHandler(
    filename = "calibration_template.csv",
    content = function(file) {
      write.csv2(data.frame(Concentration = c(0, 10, 25, 50, 100, 200, 500),
                            Response = c(0, 1000, 2500, 5000, 10000, 20000, 50000)),
                 file, row.names = FALSE)
    }
  )

  output$download_template_sample <- downloadHandler(
    filename = "sample_template.csv",
    content = function(file) {
      write.csv2(data.frame(
        Sample_Name = paste0("Sample_", 1:10),
        Response = c(1500, 3200, 4800, 2100, 5500, 1800, 4200, 3800, 2900, 4500),
        Theoretical_Concentration = c(15, 32, 48, 21, 55, 18, 42, 38, 29, 45),
        Dilution_Factor = 1,
        Class = rep(c("Control", "Treated"), each = 5),
        Batch = rep(c("A", "B"), times = 5)
      ), file, row.names = FALSE)
    }
  )

  output$download_results <- downloadHandler(
    filename = function() paste0("sample_concentrations_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- display_results()
      con <- file(file, open = "w", encoding = "UTF-8")
      on.exit(close(con))
      writeChar("\ufeff", con, eos = NULL)  # BOM so that Excel reads the accents/symbols
      if (input$csv_format == "fr") write.csv2(d, con, row.names = FALSE)
      else write.csv(d, con, row.names = FALSE)
    }
  )

  output$download_calib_plot <- downloadHandler(
    filename = function() paste0("calibration_curve_", Sys.Date(), ".pdf"),
    content = function(file) ggsave(file, plot = calib_plot_obj(), width = 10, height = 6.5,
                                    device = pdf_device())
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("concentration_plot_", Sys.Date(), ".pdf"),
    content = function(file) ggsave(file, plot = stats_plot_obj(), width = 10, height = 6.5,
                                    device = pdf_device())
  )
}

shinyApp(ui = ui, server = server)
