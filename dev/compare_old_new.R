# Compare CmbCor of the new main.R with the old Shiny operator (commit ecfbdb7).
#
#   docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/compare_old_new.R [out_dir]
#
# The old operator is replicated offline from dev/reference/:
# - pgcombat_shiny_ecfbdb7.R: the old ComBat code, unchanged
# - server_shiny_ecfbdb7.R:   the old Shiny server (reference only); its data path is
#   reproduced in old_one_step() / old_two_step() below (reshape2::acast as in server.R)
# Old one-step  = fit on all samples, output the correction computed inside fit().
# Old two-step  = step 1 fit on a view with only the fit-condition samples (Return link),
#                 step 2 apply that model to all samples (Apply saved model).
# "old + /" is the old code with only the missing "/" in fit() restored (bug D1).
#
# Output (out_dir, default dev/compare): comparison.csv (one row per case) and
# old_vs_new_scatter.png (new vs old CmbCor per case); table on stdout.

if (dir.exists("/src")) setwd("/src")
source("dev/mock_ctx.R")
if (!file.exists("dev/data/simulated.csv")) source("dev/simulate_data.R")

# reshape2 is not an operator dependency; install it into a temporary library.
if (!requireNamespace("reshape2", quietly = TRUE)) {
  lib = file.path(tempdir(), "lib"); dir.create(lib)
  install.packages("reshape2", lib = lib, repos = "https://cloud.r-project.org", quiet = TRUE)
  .libPaths(c(lib, .libPaths()))
}
library(reshape2)

old_src = readLines("dev/reference/pgcombat_shiny_ecfbdb7.R", warn = FALSE)
load_pgcombat = function(src) { e = new.env(); eval(parse(text = src), envir = e); e$pgCombat }
OldCombat = load_pgcombat(old_src)
fix_line = grep("bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,] %*% gamma.star))", old_src, fixed = TRUE)
stopifnot(length(fix_line) == 1)
slash_src = old_src
slash_src[fix_line] = paste(sub("\\s+$", "", slash_src[fix_line]), "/")
SlashCombat = load_pgcombat(slash_src)

# As server.R: df = (.y, .ri, .ci, .color = first colour factor).
old_matrices = function(df) {
  X0 = acast(df, .ri ~ .ci, value.var = ".y")
  bv = acast(df, .ri ~ .ci, value.var = ".color")[1, ]
  bv = droplevels(factor(bv))
  rowSeq = acast(df, .ri ~ .ci, value.var = ".ri")[, 1]
  colSeq = acast(df, .ri ~ .ci, value.var = ".ci")[1, ]
  dimnames(X0) = list(rowSeq = rowSeq, colSeq = colSeq)
  list(X0 = X0, bv = bv)
}
melt_result = function(Xc, X0) {
  dimnames(Xc) = dimnames(X0)
  out = melt(Xc, value.name = "CmbCor")
  data.frame(.rids = as.integer(out$rowSeq), .cids = as.integer(out$colSeq), CmbCor = out$CmbCor)
}
old_fit = function(K, m, ref, mean_only) {
  cmod = K$new()
  if (!is.null(ref)) cmod$fit(m$X0, m$bv, ref.batch = ref, mean.only = mean_only)
  else cmod$fit(m$X0, m$bv, mean.only = mean_only)
}
old_one_step = function(K, df, ref, mean_only) {
  m = old_matrices(df)
  melt_result(old_fit(K, m, ref, mean_only)$Xc, m$X0)
}
old_two_step = function(K, df, fit_ci, ref, mean_only) {
  model = old_fit(K, old_matrices(df[df$.ci %in% fit_ci, ]), ref, mean_only)
  m = old_matrices(df)
  melt_result(model$apply(m$X0, m$bv), m$X0)
}

joined = function(new, old) {
  j = merge(new, old, by = c(".rids", ".cids"), suffixes = c(".new", ".old"))
  stopifnot(nrow(j) == nrow(new), nrow(j) == nrow(old))
  j
}
max_diff = function(new, old) {
  if (is.character(old)) return(NA_real_)
  j = joined(new, old)
  max(abs(j$CmbCor.new - j$CmbCor.old))
}
try_old = function(expr) tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))

refs = list(qc = "641031403", condfit = "run01", sim = "A")
fit_cond = list(condfit = list(factors = "Grouping", values = "REF", column = "js0.Grouping"))

cases = list()
for (dsn in names(datasets)) for (model in c("L/S", "L")) for (use_ref in c(FALSE, TRUE))
  cases[[length(cases) + 1]] = list(ds = dsn, model = model, ref = use_ref, mode = "one-step")
for (model in c("L/S", "L")) for (use_ref in c(FALSE, TRUE))
  cases[[length(cases) + 1]] = list(ds = "condfit", model = model, ref = use_ref, mode = "two-step / fit condition")

out_dir = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "dev/compare"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
data_cache = lapply(datasets, read_dataset)
rows = list()
pairs = list()
for (cs in cases) {
  ds = datasets[[cs$ds]]
  props = list(ModelType = cs$model)
  ref = if (cs$ref) refs[[cs$ds]] else NULL
  if (cs$ref) props = c(props, UseReferenceBatch = "true", ReferenceBatch = ref)
  two_step = cs$mode != "one-step"
  if (two_step) props = c(props, UseFitCondition = "true", FitConditionFactors = fit_cond[[cs$ds]]$factors,
                          FitConditionValues = fit_cond[[cs$ds]]$values)
  new = run_main(ds, props, d = data_cache[[cs$ds]])
  if (!is.null(new$error)) stop("new operator failed on ", cs$ds, ": ", new$error)

  ctx = new$ctx
  df = ctx$data[, c(".y", ".ri", ".ci")]
  df$.color = ctx$data[[ctx$colors[[1]]]]
  mean_only = cs$model == "L"
  run_old = function(K) try_old(if (two_step) {
    fit_ci = which(ctx$cselect()[[fit_cond[[cs$ds]]$column]] == fit_cond[[cs$ds]]$values) - 1L
    old_two_step(K, df, fit_ci, ref, mean_only)
  } else old_one_step(K, df, ref, mean_only))
  old = run_old(OldCombat)
  slash = run_old(SlashCombat)
  case_label = sprintf("%s | %s | %s | ref: %s", cs$ds, cs$mode, cs$model, if (cs$ref) ref else "-")
  if (!is.character(old)) {
    j = joined(new$cmbcor, old)
    pairs[[length(pairs) + 1]] = data.frame(case = case_label, old = j$CmbCor.old, new = j$CmbCor.new)
  }
  rows[[length(rows) + 1]] = data.frame(
    data = cs$ds, mode = cs$mode, model = cs$model, reference = if (cs$ref) ref else "-",
    max_diff_vs_old = max_diff(new$cmbcor, old),
    max_diff_vs_old_slash = max_diff(new$cmbcor, slash),
    old_error = if (is.character(old)) old else "")
}
res = do.call(rbind, rows)
res$equal_to_old = !is.na(res$max_diff_vs_old) & res$max_diff_vs_old <= 1e-10
write.csv(res, file.path(out_dir, "comparison.csv"), row.names = FALSE)

pairs = do.call(rbind, pairs)
pairs$case = factor(pairs$case, levels = unique(pairs$case))
scatter = ggplot(pairs, aes(old, new)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 0.6, alpha = 0.4, colour = "#1f4e79") +
  facet_wrap(~case, ncol = 4, scales = "free") +
  labs(x = "CmbCor, old Shiny operator (ecfbdb7)", y = "CmbCor, new operator") +
  theme_bw(base_size = 9) + theme(strip.text = element_text(size = 7), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "old_vs_new_scatter.png"), scatter, width = 12, height = 12, dpi = 120, bg = "white")
options(width = 200)
print(res, row.names = FALSE, digits = 3)
cat("\n", sum(res$equal_to_old), "of", nrow(res), "cases equal to the old operator (max abs diff <= 1e-10)\n")
