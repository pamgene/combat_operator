# Tests of the ComBat implementation R/pgcombat.R on all test data sets.
#
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/test_pgcombat.R [out_dir]
#
# Per data set x model (L/S, L) x reference batch (none, first batch):
#   fit_vs_apply   correcting the fitted samples with apply() equals the correction inside fit()
#   ref_unchanged  reference-batch samples are returned unchanged
#   no_nan         no NaN / NA in the corrected values
#   vs_sva         max abs difference to sva::ComBat (same settings, par.prior = TRUE);
#                  needs sva in dev/.rlib (see CLAUDE.md), skipped otherwise
# plus: apply() on data with a batch the model does not know stops with a clear message.
# Exit status 1 if any check fails. Writes test_pgcombat.csv to out_dir (default dev/compare).

if (dir.exists("/src")) setwd("/src")
if (dir.exists("dev/.rlib")) .libPaths(c("dev/.rlib", .libPaths()))
source("dev/mock_ctx.R")
source("R/functions.R")
source("R/pgcombat.R")
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")
has_sva = requireNamespace("sva", quietly = TRUE)
if (!has_sva) message("sva not installed: vs_sva skipped")

tol = 1e-10       # fit vs apply, reference unchanged
tol_sva = 1e-6    # vs sva::ComBat

rows = list()
for (dsn in names(datasets)) {
  ds = datasets[[dsn]]
  ctx = mock_ctx(read_dataset(ds), ds)
  X = build_matrix(ctx$data$.ri, ctx$data$.ci, ctx$data$.y, nrow(ctx$rselect()), nrow(ctx$cselect()))
  b = factor(sample_batches(ctx$data$.ci, batch_labels(ctx$select(ctx$colors)), ncol(X)))
  for (model in c("L/S", "L")) for (ref in list(NULL, levels(b)[1])) {
    mo = model == "L"
    m = pgCombat$new()$fit(X, b, ref.batch = ref, mean.only = mo)
    ap = m$apply(X, b)
    vs_sva = NA_real_
    if (has_sva) {
      sv = suppressMessages(sva::ComBat(X, b, mod = NULL, par.prior = TRUE, mean.only = mo, ref.batch = ref))
      vs_sva = max(abs(sv - m$Xc))
    }
    rows[[length(rows) + 1]] = data.frame(
      data = dsn, model = model, reference = if (is.null(ref)) "-" else ref,
      fit_vs_apply = max(abs(m$Xc - ap)),
      ref_unchanged = if (is.null(ref)) NA else max(abs(m$Xc[, b == ref] - X[, b == ref])) <= tol,
      no_nan = !anyNA(m$Xc) && !anyNA(ap),
      vs_sva = vs_sva)
  }
}
res = do.call(rbind, rows)
res$pass = res$fit_vs_apply <= tol & (is.na(res$ref_unchanged) | res$ref_unchanged) & res$no_nan &
  (is.na(res$vs_sva) | res$vs_sva <= tol_sva)

# apply() with a batch the model does not know
ds = datasets$sim; ctx = mock_ctx(read_dataset(ds), ds)
X = build_matrix(ctx$data$.ri, ctx$data$.ci, ctx$data$.y, nrow(ctx$rselect()), nrow(ctx$cselect()))
b = factor(sample_batches(ctx$data$.ci, batch_labels(ctx$select(ctx$colors)), ncol(X)))
m = pgCombat$new()$fit(X[, b != "C"], droplevels(b[b != "C"]))
msg = tryCatch({ m$apply(X, b); "no error" }, error = function(e) conditionMessage(e))
unknown_ok = grepl("C", msg) && grepl("batch", msg, ignore.case = TRUE) && msg != "no error"

options(width = 200)
print(res, row.names = FALSE, digits = 3)
cat("\napply() with unknown batch C:", msg, "\n")
cat("unknown-batch check:", if (unknown_ok) "pass" else "FAIL", "\n")
out_dir = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/compare"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "test_pgcombat.csv"), row.names = FALSE)
n_fail = sum(!res$pass) + !unknown_ok
cat(if (n_fail == 0) "ALL PASS" else paste(n_fail, "FAILED"), "\n")
quit(status = as.integer(n_fail > 0))
