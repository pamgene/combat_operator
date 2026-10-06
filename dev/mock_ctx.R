# Mock of the Tercen ctx API used by main.R, plus the three test data sets with
# their crosstab projection. Sourced by dev/preview.R and dev/compare_old_new.R.
#
# Tercen indexes row/column factor combinations in sorted order, starting at 0.

suppressMessages({library(dplyr); library(ggplot2); library(R6)})

# Data sets: file, crosstab projection (factor names as in the CSV header).
datasets = list(
  qc = list(file = "dev/data/example_data_QC.csv",
            rows = "S100QC..ID",
            cols = c("S100QC..Barcode", "S100QC..Row", "S100QC..Supergroup", "S100QC..Test Condition"),
            colours = "S100QC..Barcode",
            y = "S100QC.logTransformed"),
  condfit = list(file = "dev/data/example_data_conditionfit_ref.csv",
                 rows = "VSN-QC..ID",
                 cols = c("VSN-QC..Barcode", "VSN-QC..Array", "js0.Grouping", "js0.Sample.name"),
                 colours = "js0.Run",
                 y = "VSN-QC.identity"),
  sim = list(file = "dev/data/simulated.csv",
             rows = "sim.peptide",
             cols = c("sim.sample", "sim.condition"),
             colours = "sim.batch",
             y = "sim.value")
)

read_dataset = function(ds) {
  d = read.csv(ds$file, check.names = FALSE, stringsAsFactors = FALSE,
               colClasses = "character")
  d[[ds$y]] = as.numeric(d[[ds$y]])
  d
}

# Every property gets its operator.json defaultValue, as Tercen passes it.
json_defaults = local({
  props = jsonlite::fromJSON("operator.json")$properties
  setNames(as.list(as.character(props$defaultValue)), props$name)
})

mock_ctx = function(d, ds, props = list(), colours = ds$colours) {
  props = modifyList(json_defaults, props)
  col_df = d[, ds$cols, drop = FALSE] %>% distinct() %>% arrange(across(everything()))
  row_df = d[, ds$rows, drop = FALSE] %>% distinct() %>% arrange(across(everything()))
  col_key = do.call(paste, c(col_df, sep = "\r"))
  row_key = do.call(paste, c(row_df, sep = "\r"))
  qt = data.frame(
    .ri = match(do.call(paste, c(d[, ds$rows, drop = FALSE], sep = "\r")), row_key) - 1L,
    .ci = match(do.call(paste, c(d[, ds$cols, drop = FALSE], sep = "\r")), col_key) - 1L,
    .y = d[[ds$y]]
  )
  qt = cbind(qt, d[, colours, drop = FALSE])
  ctx = list(
    colors = colours,
    cnames = ds$cols,
    rnames = ds$rows,
    op.value = function(name, type, default) if (!is.null(props[[name]])) type(props[[name]]) else default,
    select = function(cols) qt[, cols, drop = FALSE],
    cselect = function() col_df,
    rselect = function() row_df,
    addNamespace = function(x) x,
    data = qt
  )
  class(ctx) = "mockctx"
  ctx
}
select.mockctx = function(.data, ...) dplyr::select(.data$data, ...)

# Intercept main.R's output chain: keep the CmbCor table and the PNG path.
as_relation = function(x, ...) x
left_join_relation = function(x, ...) x
as_join_operator = function(x, ...) x
file_to_tercen = function(file) file
save_relation = function(x, ctx) assign("captured", list(cmbcor = x[[1]], png = x[[2]]), envir = globalenv())

main_src = local({
  e = parse("main.R")
  drop = vapply(e, function(x) identical(x, quote(library(tercen))), logical(1))
  e[!drop]
})

# Run the unmodified main.R on a data set with the given properties.
# Returns list(cmbcor = data.frame(.rids, .cids, CmbCor), png, ctx) or list(error = message).
run_main = function(ds, props = list(), colours = ds$colours, d = read_dataset(ds)) {
  ctx_obj = mock_ctx(d, ds, props, colours)
  assign("tercenCtx", function() ctx_obj, envir = globalenv())
  rm(list = intersect("captured", ls(globalenv())), envir = globalenv())
  err = tryCatch({ eval(main_src, envir = globalenv()); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(err)) return(list(error = err))
  c(get("captured", envir = globalenv()), list(ctx = ctx_obj))
}
