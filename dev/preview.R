# Preview of main.R on the test data sets, without Tercen.
#
# Runs the unmodified main.R against a mock Tercen context (dev/mock_ctx.R).
# Run it inside the operator image so R and packages match production:
#
#   docker build -t combat_operator:dev .
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/preview.R
#
# Output: dev/preview/<case>.png, dev/preview/<case>.csv (CmbCor), dev/preview/errors.txt,
# dev/preview/index.html (overview).

if (dir.exists("/src")) setwd("/src")
source("dev/mock_ctx.R")
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")

out_dir = "dev/preview"
dir.create(out_dir, showWarnings = FALSE)
unlink(list.files(out_dir, full.names = TRUE))

cond_ref = list(UseFitCondition = "true", FitConditionFactors = "Grouping", FitConditionValues = "REF")
sim_ctl = list(UseFitCondition = "true", FitConditionFactors = "condition", FitConditionValues = "Control")

# name, data set, properties, colour factors (NULL = data set default)
cases = list(
  list("qc_LS", "qc", list()),
  list("qc_L", "qc", list(ModelType = "L")),
  list("qc_LS_ref", "qc", list(UseReferenceBatch = "true", ReferenceBatch = "641031403")),
  list("condfit_LS", "condfit", list()),
  list("condfit_LS_fitREF", "condfit", cond_ref),
  list("condfit_L_fitREF", "condfit", c(cond_ref, ModelType = "L")),
  list("condfit_LS_fitREF_ref", "condfit", c(cond_ref, UseReferenceBatch = "true", ReferenceBatch = "run01")),
  list("qc_L_2colours_ref", "qc", list(ModelType = "L", UseReferenceBatch = "true", ReferenceBatch = "641031404 ; Sgroup1"),
       c("S100QC..Barcode", "S100QC..Supergroup")),
  list("sim_LS", "sim", list()),
  list("sim_L", "sim", list(ModelType = "L")),
  list("sim_LS_ref", "sim", list(UseReferenceBatch = "true", ReferenceBatch = "A")),
  list("sim_LS_fitControl", "sim", sim_ctl),
  list("sim_LS_fitControl_ref", "sim", c(sim_ctl, UseReferenceBatch = "true", ReferenceBatch = "A")),
  # error cases: the message is what the user sees
  list("ERR_qc_fitControl_LS", "qc", list(UseFitCondition = "true", FitConditionFactors = "Test Condition", FitConditionValues = "Control")),
  list("ERR_qc_fitControl_L", "qc", list(ModelType = "L", UseFitCondition = "true", FitConditionFactors = "Test Condition", FitConditionValues = "Control")),
  list("ERR_unknown_factor", "condfit", list(UseFitCondition = "true", FitConditionFactors = "Supergroup", FitConditionValues = "REF")),
  list("ERR_value_count", "condfit", list(UseFitCondition = "true", FitConditionFactors = "Grouping", FitConditionValues = "REF;X")),
  list("ERR_no_match", "condfit", list(UseFitCondition = "true", FitConditionFactors = "Grouping", FitConditionValues = "CTRL")),
  list("ERR_ref_unknown", "qc", list(UseReferenceBatch = "true", ReferenceBatch = "123")),
  list("ERR_2colours_fitREF_too_fine", "condfit", c(cond_ref, ModelType = "L"), c("js0.Run", "VSN-QC..Barcode")),
  list("ERR_ref_2colours", "condfit", list(UseReferenceBatch = "true", ReferenceBatch = "run01"), c("js0.Run", "VSN-QC..Barcode"))
)

data_cache = lapply(datasets, read_dataset)
errors = character(0)
for (cs in cases) {
  name = cs[[1]]; ds = datasets[[cs[[2]]]]
  colours = if (length(cs) > 3) cs[[4]] else ds$colours
  res = run_main(ds, cs[[3]], colours, data_cache[[cs[[2]]]])
  if (!is.null(res$error)) {
    errors = c(errors, sprintf("%s:\n  %s\n", name, res$error))
    cat(name, ": ERROR\n")
    next
  }
  file.copy(res$png, file.path(out_dir, paste0(name, ".png")), overwrite = TRUE)
  write.csv(res$cmbcor, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat(name, ": ok, NaN =", sum(is.nan(res$cmbcor$CmbCor)), "\n")
}
writeLines(errors, file.path(out_dir, "errors.txt"))

pngs = sort(basename(list.files(out_dir, pattern = "[.]png$")))
writeLines(c('<!doctype html><meta charset="utf-8"><title>ComBat preview</title>',
             '<style>body{font-family:sans-serif;margin:16px;background:#fff}',
             'figure{margin:0 0 32px}figcaption{font-weight:bold;margin-bottom:6px}',
             'img{max-width:100%;border:1px solid #ddd}pre{white-space:pre-wrap}</style>',
             '<h2>Error messages</h2><pre>', gsub("<", "&lt;", errors), '</pre>',
             sprintf('<figure><figcaption>%s</figcaption><img src="%s"></figure>',
                     sub("[.]png$", "", pngs), pngs)),
           file.path(out_dir, "index.html"))
cat("wrote", length(pngs), "PNGs and", length(errors), "error messages to", out_dir, "\n")
