# ============================================================================
# MRDS duplicate matcher: Shiny app
# ----------------------------------------------------------------------------
# Identifies duplicate detections between two tandem observers (A = front,
# B = rear) on an aerial survey and builds the 11 / 10 / 01 detection history
# needed for mark-recapture distance sampling.
#
# Runs as an ordinary Shiny app (shiny::runApp()) and, exported with
# shinylive, entirely in the browser as WebAssembly. Functions are in R/.
# ============================================================================

library(shiny)
library(bslib)
library(DT)
# Loaded by R/01_mrds_duplicate_match.R as well; listed here so shinylive
# bundles them.
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(geosphere)
library(clue)

# Shiny sources R/ automatically; this covers source("app.R") and tests.
if (!exists("run_matching", mode = "function"))
  for (f in sort(list.files("R", full.names = TRUE, pattern = "\\.R$"))) source(f)

options(shiny.maxRequestSize = 50 * 1024^2)

METHODS <- c("Fellegi-Sunter (EM)" = "fellegi_sunter",
             "Discrepancy mixture (Hamilton et al.)" = "discrepancy",
             "Fixed tolerance" = "fixed")
ASSIGNMENTS <- c("Hungarian, one-to-one" = "hungarian",
                 "Stochastic clusters (Hamilton et al.)" = "stochastic_cluster",
                 "Metropolis-Hastings sample" = "mh_sample")
SPACE_FIELDS <- c("Automatic" = "auto", "Sighting angle" = "angle",
                  "Perpendicular distance" = "distance", "GPS separation (m)" = "d_m")

hint <- function(...) tags$p(class = "text-muted small mb-2", ...)

# A tab whose contents keep their natural height (a plain div stops bslib
# stretching and squashing tables to fit the card).
panel <- function(title, ...) nav_panel(title, div(...))

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

sidebar_ui <- sidebar(
  width = 370,
  accordion(
    open = c("Data", "Scoring"), multiple = TRUE,

    accordion_panel(
      "Data",
      radioButtons("source", NULL, inline = TRUE,
                   c("Simulated example" = "example", "Upload CSV" = "upload")),
      conditionalPanel(
        "input.source == 'example'",
        hint("A simulated six-transect survey with known duplicates, so precision",
             "and recall can be shown."),
        downloadLink("dl_example", "Download the example as CSV")),
      conditionalPanel(
        "input.source == 'upload'",
        fileInput("file", NULL, accept = c(".csv", "text/csv"),
                  placeholder = "One row per detection"),
        hint("The file is read in your browser and is not sent anywhere."),
        uiOutput("mapping_ui"))
    ),

    accordion_panel(
      "Gate",
      hint("Only A-B pairs inside the gate are compared. Set it several times wider",
           "than the disagreement you expect between true duplicates."),
      numericInput("gate_t", "Maximum time gap (s)", 15, min = 0.1, step = 1),
      numericInput("gate_s", "Maximum separation (m)", 300, min = 1, step = 10)
    ),

    accordion_panel(
      "Scoring",
      selectInput("method", "Method", METHODS),
      conditionalPanel(
        "input.method == 'fellegi_sunter'",
        hint("Starting values for the EM algorithm."),
        numericInput("init_sigma_t", "SD of time offset (s)", 3, min = 0.1, step = 0.5),
        numericInput("init_sigma_s", "Position scale (m)", 50, min = 1, step = 5),
        numericInput("init_sigma_g", "SD of log group-size ratio", 0.3, min = 0.01, step = 0.05),
        numericInput("max_iter", "Maximum iterations", 200, min = 1, step = 50)),
      conditionalPanel(
        "input.method == 'discrepancy'",
        hint("Discrepancy = min(|dt|, cap) / time scale + spatial difference / spatial",
             "scale. The scales need tuning to the survey."),
        numericInput("time_cap", "Time cap (s)", 25, min = 0.1),
        numericInput("time_scale", "Time scale (s)", 25, min = 0.1),
        selectInput("disc_space_field", "Spatial field", SPACE_FIELDS),
        numericInput("space_scale", "Spatial scale (blank = default)", NA, min = 0),
        hint("Default: 65 for angles, otherwise the separation gate."),
        numericInput("discrep_max", "Maximum discrepancy", 1, min = 0.01, step = 0.1),
        checkboxInput("centre_time", "Remove the median time lag first", FALSE),
        selectInput("known_nonmatch", "Known non-matches",
                    c("None" = "none", "Same-observer pairs" = "same_observer",
                      "Opposite-side pairs" = "opposite_side"))),
      conditionalPanel(
        "input.method == 'fixed'",
        numericInput("time_tol", "Time tolerance (s)", 10, min = 0),
        numericInput("time_offset", "Time offset, B minus A (s)", 0, step = 0.5),
        selectInput("fix_space_field", "Spatial field", SPACE_FIELDS),
        numericInput("space_tol", "Spatial tolerance", 15, min = 0),
        hint("In the units of the spatial field: degrees for angles, your distance",
             "units, or metres for GPS separation."))
    ),

    accordion_panel(
      "Assignment",
      selectInput("assignment", "Method", ASSIGNMENTS),
      conditionalPanel(
        "input.assignment != 'stochastic_cluster'",
        sliderInput("threshold", "Match probability threshold", 0.05, 0.95, 0.5, 0.05)),
      conditionalPanel(
        "input.assignment != 'hungarian'",
        numericInput("K", "Number of draws", 20, min = 1, max = 500, step = 1)),
      conditionalPanel(
        "input.assignment == 'mh_sample'",
        numericInput("n_burn", "Burn-in iterations", 2000, min = 0, step = 500),
        numericInput("n_thin", "Thinning interval", 50, min = 1, step = 10),
        hint("Per transect when a transect column is mapped. Each iteration proposes",
             "a change to one detection, so allow several per detection.")),
      numericInput("seed", "Random seed (blank = none)", 1, step = 1)
    )
  ),
  actionButton("run", "Run matching", class = "btn-primary w-100 mt-2")
)

ui <- page_sidebar(
  title = "MRDS duplicate matcher",
  theme = bs_theme(version = 5, primary = "#1f5a7a"),
  sidebar = sidebar_ui,
  fillable = FALSE,   # let the page scroll so tables keep their natural height
  navset_card_tab(
    id = "tabs",
    panel(
      "Data",
      uiOutput("data_notes"),
      DTOutput("data_tbl")),
    panel(
      "Summary",
      uiOutput("summary_boxes"),
      uiOutput("fit_status"),
      layout_columns(
        col_widths = c(6, 6),
        div(h6("Estimated parameters"), tableOutput("param_tbl")),
        div(h6("Detection histories"), tableOutput("hist_counts"),
            uiOutput("truth_eval")))),
    panel(
      "Diagnostics",
      layout_columns(
        col_widths = c(6, 6),
        plotOutput("plot_p", height = 320),
        plotOutput("plot_scatter", height = 320)),
      conditionalPanel("output.is_discrepancy", plotOutput("plot_discrep", height = 300)),
      layout_columns(
        col_widths = c(5, 7),
        div(h6("Candidate pairs by match probability"), tableOutput("class_tbl")),
        div(h6("Residual correlations within the match class"),
            tableOutput("dep_tbl"),
            hint("Values near zero support the conditional-independence assumption.",
                 "dt = time offset, d_m = separation, lsr = log size ratio, sp = species.")))),
    panel(
      "Matches",
      uiOutput("draw_ui"),
      div(class = "mb-2",
          downloadButton("dl_matches", "Matched pairs (CSV)", class = "btn-sm"),
          downloadButton("dl_draws", "All draws (CSV)", class = "btn-sm"),
          downloadButton("dl_pairs", "All candidate pairs (CSV)", class = "btn-sm")),
      DTOutput("match_tbl"),
      uiOutput("freq_ui")),
    panel(
      "Detection history",
      hint("One row per group: 11 = seen by both, 10 = A only, 01 = B only.",
           "Follows the draw selected on the Matches tab."),
      div(class = "mb-2",
          downloadButton("dl_history", "Detection history (CSV)", class = "btn-sm")),
      DTOutput("history_tbl")),
    panel(
      "About",
      div(style = "max-width: 760px",
        h5("What this does"),
        p("Two observers on the same side of the aircraft record detections independently.",
          "Before a mark-recapture distance sampling (MRDS) detection function can be",
          "fitted, records from A and B that refer to the same group have to be identified.",
          "This app scores every plausible A-B pair, assigns duplicates, and writes the",
          "detection history."),
        h5("Input"),
        p("A CSV with one row per detection. Required: observer, time, longitude and",
          "latitude in decimal degrees (WGS84), and group size. Optional: a detection ID,",
          "transect, side of aircraft, species, perpendicular distance and sighting angle.",
          "Time can be numeric seconds, a clock time (HH:MM:SS) or a date-time. Pairs are",
          "only formed within a transect and side when those columns are mapped."),
        h5("Scoring methods"),
        tags$ul(
          tags$li(strong("Fellegi-Sunter (EM)."), "A two-class record-linkage model fitted by EM.",
                  "For true duplicates the time offset is Normal, the separation Rayleigh, the",
                  "log group-size ratio Normal about zero and species agreement Bernoulli;",
                  "non-matches follow what random pairing inside the gate would give."),
          tags$li(strong("Discrepancy mixture."), "Hamilton et al. (2018): a scaled sum of the time",
                  "and spatial differences, modelled as a two-component truncated-Normal mixture."),
          tags$li(strong("Fixed tolerance."), "A pair is a duplicate if both differences are",
                  "inside fixed tolerances.")),
        h5("Assignment methods"),
        tags$ul(
          tags$li(strong("Hungarian."), "The single most probable one-to-one matching, keeping",
                  "pairs above the probability threshold."),
          tags$li(strong("Stochastic clusters."), "Hamilton et al.: each pair is selected with its",
                  "match probability and clusters are resolved by anchoring and merging, so one",
                  "observer's sighting can match a group the other observer split in two."),
          tags$li(strong("Metropolis-Hastings sample."), "Draws one-to-one matchings from their",
                  "posterior. Use the draws for multiple imputation in the MRDS analysis.")),
        h5("Limits"),
        p("Everything runs in your browser. It is slower than desktop R, and very large",
          "files may run out of memory; map a transect column to keep the comparisons",
          "manageable. Abundance estimation is not part of this app."),
        h5("Reference"),
        p("Hamilton, O. N. P., Kincaid, S. E., Constantine, R., Kozmian-Ledward, L., Walker,",
          "C. G. & Fewster, R. M. (2018) Accounting for uncertainty in duplicate identification",
          "and group size judgments in mark-recapture distance sampling.",
          em("Methods in Ecology and Evolution"), "9: 354-362.")))
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  example_data <- simulate_example()

  # ---- data ---------------------------------------------------------------
  raw <- reactive({
    req(input$source == "upload", input$file)
    d <- tryCatch(
      read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE,
               na.strings = c("", "NA", "na", "N/A"), strip.white = TRUE),
      error = function(e) NULL)
    validate(need(!is.null(d) && nrow(d) > 0 && ncol(d) > 1,
                  "Could not read that file as a CSV with a header row."))
    names(d) <- make.unique(trimws(names(d)))
    d
  })

  output$mapping_ui <- renderUI({
    d <- raw()
    cols <- names(d)
    sel <- function(role) {
      r <- COLUMN_ROLES[[role]]
      selectInput(paste0("map_", role),
                  if (r$required) r$label else paste(r$label, "(optional)"),
                  choices = if (r$required) c("Choose a column" = "", cols) else c(NONE, cols),
                  selected = { g <- guess_column(role, cols)
                               if (g == NONE && r$required) "" else g })
    }
    tagList(
      h6("Columns"),
      sel("observer"),
      uiOutput("observer_values_ui"),
      sel("time"),
      selectInput("time_format", "Time format",
                  c("Detect automatically" = "auto", "Seconds (numeric)" = "seconds",
                    "Clock time (HH:MM:SS)" = "clock", "Date and time" = "datetime")),
      sel("lon"), sel("lat"), sel("group_size"),
      sel("id"), sel("transect"), sel("side"), sel("species"),
      sel("distance"), sel("angle"))
  })

  output$observer_values_ui <- renderUI({
    d <- raw(); req(input$map_observer, input$map_observer %in% names(d))
    vals <- sort(unique(as.character(d[[input$map_observer]])))
    vals <- vals[!is.na(vals)]
    validate(need(length(vals) >= 2, "The observer column needs at least two values."),
             need(length(vals) <= 50, "That column has too many values to be the observer."))
    tagList(
      selectInput("obs_A", "Observer A (front)", vals, vals[1]),
      selectInput("obs_B", "Observer B (rear)", vals, vals[2]))
  })

  std <- reactive({
    if (input$source == "example")
      return(list(data = example_data,
                  notes = "Simulated data. The truth_id column gives the true group."))
    d <- raw()
    req(input$map_observer, input$obs_A, input$obs_B)
    map <- lapply(setNames(names(COLUMN_ROLES), names(COLUMN_ROLES)),
                  function(r) input[[paste0("map_", r)]])
    out <- tryCatch(standardise_detections(d, map, input$obs_A, input$obs_B,
                                           input$time_format %||% "auto"),
                    error = function(e) conditionMessage(e))
    validate(need(is.list(out), if (is.character(out)) out else "Could not read the data."))
    out
  })

  output$data_notes <- renderUI({
    s <- std(); d <- s$data
    n_tr <- if ("transect" %in% names(d)) length(unique(d$transect)) else NA
    tagList(
      p(strong(sprintf("%d detections: %d by A, %d by B%s.", nrow(d),
                       sum(d$observer == "A"), sum(d$observer == "B"),
                       if (is.na(n_tr)) "" else sprintf(", on %d transects", n_tr)))),
      if (length(s$notes)) tags$ul(class = "small text-muted", lapply(s$notes, tags$li)))
  })

  output$data_tbl <- renderDT({
    s <- tryCatch(std(), error = function(e) NULL)   # message shown once, above
    req(s)
    datatable(s$data, rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })

  output$dl_example <- downloadHandler(
    filename = "example_tandem_survey.csv",
    content = function(file) write.csv(example_data, file, row.names = FALSE))

  # ---- run ----------------------------------------------------------------
  num_or <- function(x, default) if (is.null(x) || is.na(x)) default else x

  result <- eventReactive(input$run, {
    det <- std()$data
    method_args <- switch(input$method,
      fellegi_sunter = list(
        init = list(mu_t = 0, sigma_t = num_or(input$init_sigma_t, 3),
                    sigma_s = num_or(input$init_sigma_s, 50),
                    sigma_g = num_or(input$init_sigma_g, 0.3),
                    eps = 0.02, pi_match = NULL),
        max_iter = num_or(input$max_iter, 200)),
      discrepancy = list(
        time_cap = num_or(input$time_cap, 25), time_scale = num_or(input$time_scale, 25),
        space_scale = if (is.na(num_or(input$space_scale, NA))) NULL else input$space_scale,
        space_field = input$disc_space_field,
        discrep_max = num_or(input$discrep_max, 1),
        centre_time = isTRUE(input$centre_time)),
      fixed = list(
        time_tol = num_or(input$time_tol, 10), space_tol = num_or(input$space_tol, 15),
        space_field = input$fix_space_field, time_offset = num_or(input$time_offset, 0)))

    withProgress(message = "Matching", value = 0.1, {
      tryCatch(
        run_matching(det, method = input$method, assignment = input$assignment,
                     max_time_gap_s = num_or(input$gate_t, 15),
                     max_space_gap_m = num_or(input$gate_s, 300),
                     method_args = method_args,
                     known_nonmatch = if (input$method == "discrepancy") input$known_nonmatch
                                      else "none",
                     threshold = input$threshold,
                     K = if (input$assignment == "hungarian") 1 else
                           max(1, min(500, round(num_or(input$K, 20)))),
                     n_burn = max(0, round(num_or(input$n_burn, 2000))),
                     n_thin = max(1, round(num_or(input$n_thin, 50))),
                     seed = input$seed,
                     progress = function(msg) incProgress(0.25, detail = msg)),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = NULL)
          NULL
        })
    })
  })

  observeEvent(result(), {
    if (!is.null(result())) nav_select("tabs", "Summary")
  })

  res <- reactive({ r <- result(); req(r); r })

  draw_k <- reactive({
    K <- length(res()$draws)
    k <- input$draw
    if (is.null(k) || is.na(k)) 1 else max(1, min(K, round(k)))
  })

  matches  <- reactive(match_table(res(), draw_k()))
  history  <- reactive(build_history(res()$detections, res()$draws[[draw_k()]]))

  # ---- summary ------------------------------------------------------------
  output$summary_boxes <- renderUI({
    r <- res(); h <- history()
    box <- function(title, value) value_box(title, value, theme = "light")
    layout_columns(
      fill = FALSE,
      box("Candidate pairs", format(nrow(r$fit$pairs), big.mark = ",")),
      box("Duplicates", format(sum(h$seen == "11"), big.mark = ",")),
      box("A only", format(sum(h$seen == "10"), big.mark = ",")),
      box("B only", format(sum(h$seen == "01"), big.mark = ",")),
      box("Groups", format(nrow(h), big.mark = ",")))
  })

  output$fit_status <- renderUI({
    r <- res(); f <- r$fit
    msgs <- list()
    if (r$method == "fellegi_sunter") {
      msgs <- c(msgs, list(if (isTRUE(f$converged))
        div(class = "alert alert-success py-2", sprintf("EM converged in %d iterations.", f$n_iter))
      else
        div(class = "alert alert-warning py-2",
            sprintf("EM did not converge in %d iterations. Raise the maximum or change the starting values.",
                    f$n_iter))))
      msgs <- c(msgs, lapply(f$warnings, function(w) div(class = "alert alert-warning py-2", w)))
    } else if (r$method == "discrepancy") {
      msgs <- c(msgs, list(if (isTRUE(f$converged))
        div(class = "alert alert-success py-2", "Mixture model converged.")
      else
        div(class = "alert alert-warning py-2",
            paste("Mixture optimiser did not converge:", f$optim_message))))
    }
    if (length(r$draws) > 1)
      msgs <- c(msgs, list(div(class = "alert alert-info py-2",
        sprintf("Showing draw %d of %d. Counts differ between draws.", draw_k(), length(r$draws)))))
    tagList(msgs)
  })

  output$param_tbl <- renderTable(param_table(res()$fit), striped = TRUE, spacing = "s")

  output$hist_counts <- renderTable({
    h <- history()
    tibble(History = c("11 (both)", "10 (A only)", "01 (B only)"),
           Groups = c(sum(h$seen == "11"), sum(h$seen == "10"), sum(h$seen == "01")))
  }, striped = TRUE, spacing = "s")

  output$truth_eval <- renderUI({
    ev <- evaluate_truth(res()$detections, res()$draws[[draw_k()]])
    req(ev)
    tagList(
      h6("Against the simulated truth"),
      tags$table(class = "table table-sm table-striped", style = "width:auto",
        tags$tr(tags$td("True duplicates"), tags$td(ev$n_true)),
        tags$tr(tags$td("Declared duplicates"), tags$td(ev$n_declared)),
        tags$tr(tags$td("Correct"), tags$td(ev$true_pos)),
        tags$tr(tags$td("Precision"), tags$td(sprintf("%.3f", ev$precision))),
        tags$tr(tags$td("Recall"), tags$td(sprintf("%.3f", ev$recall)))))
  })

  # ---- diagnostics --------------------------------------------------------
  output$is_discrepancy <- reactive(isTRUE(result()$method == "discrepancy"))
  outputOptions(output, "is_discrepancy", suspendWhenHidden = FALSE)

  p_colours <- function(p) {
    ramp <- colorRamp(c("#9aa5ad", "#e9b44c", "#b3261e"))
    rgb(ramp(pmin(pmax(p, 0), 1)), maxColorValue = 255, alpha = 170)
  }

  output$plot_p <- renderPlot({
    p <- res()$fit$pairs$p_match
    hist(p, breaks = seq(0, 1, 0.05), col = "#1f5a7a", border = "white",
         main = "Match probability of candidate pairs",
         xlab = "P(match)", ylab = "Candidate pairs")
    if (res()$assignment != "stochastic_cluster")
      abline(v = res()$settings$threshold, lty = 2, col = "#b3261e")
  })

  output$plot_scatter <- renderPlot({
    r <- res(); pr <- r$fit$pairs
    o <- order(pr$p_match)
    plot(pr$dt[o], pr$d_m[o], pch = 16, cex = 0.8, col = p_colours(pr$p_match[o]),
         xlab = "Time offset, B minus A (s)", ylab = "Separation (m)",
         main = "Candidate pairs inside the gate",
         xlim = c(-1, 1) * r$settings$max_time_gap_s, ylim = c(0, r$settings$max_space_gap_m))
    m <- matches()
    points(m$dt, m$d_m, cex = 1.3, col = "#1f5a7a")
    legend("topleft", bty = "n", cex = 0.85, pch = c(16, 16, 16, 1),
           col = c(p_colours(c(0, 0.5, 1)), "#1f5a7a"),
           legend = c("P = 0", "P = 0.5", "P = 1", "Declared duplicate"))
  })

  output$plot_discrep <- renderPlot({
    r <- res(); req(r$method == "discrepancy")
    f <- r$fit; th <- f$params; dm <- f$settings$discrep_max
    x <- f$pairs$discrep[f$pairs$discrep <= dm]
    hist(x, breaks = 30, freq = FALSE, col = "#d6dde2", border = "white",
         main = "Discrepancy with the fitted mixture", xlab = "Discrepancy")
    tn <- function(x, m, s) dnorm(x, m, s) / (pnorm(dm, m, s) - pnorm(0, m, s))
    g <- seq(0, dm, length.out = 300)
    lines(g, th$dupP * tn(g, th$mu1, th$sig1), col = "#b3261e", lwd = 2)
    lines(g, (1 - th$dupP) * tn(g, th$mu2, th$sig2), col = "#1f5a7a", lwd = 2)
    legend("top", bty = "n", lwd = 2, col = c("#b3261e", "#1f5a7a"),
           legend = c("Near (duplicates)", "Far (non-duplicates)"))
  })

  output$class_tbl <- renderTable({
    tb <- diag_classification_table(res()$fit)
    tibble(`P(match)` = names(tb), Pairs = as.integer(tb))
  }, striped = TRUE, spacing = "s")

  output$dep_tbl <- renderTable({
    m <- tryCatch(diag_residual_dependence(res()$fit)$match_class, error = function(e) NULL)
    validate(need(!is.null(m), "Not available for this fit."))
    as.data.frame(round(m, 3)) |> tibble::rownames_to_column(" ")
  }, striped = TRUE, spacing = "s", na = "-")

  # ---- matches ------------------------------------------------------------
  output$draw_ui <- renderUI({
    K <- length(res()$draws)
    if (K > 1) numericInput("draw", sprintf("Draw (1 to %d)", K), 1, min = 1, max = K, step = 1)
  })

  round_cols <- function(dt, d, digits = 2, cols = c("dt", "d_m", "distance", "distance_A",
                                                     "distance_B")) {
    cols <- intersect(cols, names(d))
    if (length(cols)) dt <- formatRound(dt, cols, digits)
    if ("p_match" %in% names(d)) dt <- formatRound(dt, "p_match", 3)
    dt
  }

  output$match_tbl <- renderDT({
    m <- matches()
    datatable(m, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      round_cols(m)
  })

  output$freq_ui <- renderUI({
    req(length(res()$draws) > 1)
    tagList(h6(class = "mt-3", "How often each pair was declared a duplicate"),
            DTOutput("freq_tbl"))
  })

  output$freq_tbl <- renderDT({
    f <- pair_frequency(res())
    datatable(f, rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE, dom = "ftip")) |>
      formatRound(c("proportion", "p_match"), 3)
  })

  output$history_tbl <- renderDT({
    h <- history()
    lead <- intersect(c("group_id", "seen", "id_A", "id_B", "transect", "size",
                        "group_size_A", "group_size_B", "distance"), names(h))
    h <- h[, c(lead, setdiff(names(h), lead))]
    datatable(h, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
      round_cols(h)
  })

  csv_download <- function(name, data) downloadHandler(
    filename = function() sprintf("%s_%s.csv", name, format(Sys.Date(), "%Y%m%d")),
    content = function(file) write.csv(data(), file, row.names = FALSE, na = ""))

  output$dl_matches <- csv_download("matched_pairs", matches)
  output$dl_history <- csv_download("detection_history", history)
  output$dl_draws   <- csv_download("matching_draws", reactive(all_draws_long(res())))
  output$dl_pairs   <- csv_download("candidate_pairs", reactive(res()$fit$pairs))
}

shinyApp(ui, server)
