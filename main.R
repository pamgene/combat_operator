library(tercen)
library(dplyr)
library(ggplot2)
library(R6)

source("R/pgcombat.R")
source("R/functions.R")

ctx = tercenCtx()

mean_only = ctx$op.value("ModelType", as.character, "L/S") == "L"
use_ref = ctx$op.value("UseReferenceBatch", as.logical, FALSE)
ref_string = ctx$op.value("ReferenceBatch", as.character, "")
use_fit_condition = ctx$op.value("UseFitCondition", as.logical, FALSE)
fit_factors = ctx$op.value("FitConditionFactors", as.character, "")
fit_values = ctx$op.value("FitConditionValues", as.character, "")

if (length(ctx$colors) == 0) {
  stop("No batch factor: put the technical batch factor(s), e.g. Run or Barcode, on colour.", call. = FALSE)
}

df = ctx %>% select(.ri, .ci, .y)
col_df = ctx$cselect()
n_rows = nrow(ctx$rselect())
n_cols = nrow(col_df)

X0 = build_matrix(df$.ri, df$.ci, df$.y, n_rows, n_cols)
batch = sample_batches(df$.ci, batch_labels(ctx$select(ctx$colors)), n_cols)

fit_sample = if (use_fit_condition) {
  fit_condition_samples(fit_factors, fit_values, as.data.frame(col_df))
} else {
  rep(TRUE, n_cols)
}
check_fit_counts(batch[fit_sample], batch, mean_only, use_fit_condition)
ref = if (use_ref) resolve_reference(ref_string, batch, ctx$colors) else NULL

model = pgCombat$new()$fit(X0[, fit_sample, drop = FALSE], factor(batch[fit_sample]),
                           ref.batch = ref, mean.only = mean_only)
Xc = model$apply(X0, batch)
dimnames(Xc) = dimnames(X0)

params = parameter_table(mean_only, ctx$colors, ref, if (use_fit_condition) fit_factors, fit_values)
param_file = file.path(tempdir(), "combat_parameters.csv")
write.csv(params, param_file, row.names = FALSE)

plot_file = tempfile(fileext = ".png")
ggsave(plot_file, pca_plot(X0, Xc, batch, if (use_fit_condition) fit_sample, parameter_text(params)),
       width = 10, height = 5.3, dpi = 150, bg = "white")

result = data.frame(.rids = rep(0:(n_rows - 1), n_cols),
                    .cids = rep(0:(n_cols - 1), each = n_rows),
                    CmbCor = as.vector(Xc))

main_rel = result %>%
  ctx$addNamespace() %>%
  as_tibble() %>%
  as_relation() %>%
  left_join_relation(ctx$crelation, ".cids", ctx$crelation$rids) %>%
  left_join_relation(ctx$rrelation, ".rids", ctx$rrelation$rids) %>%
  as_join_operator(c(ctx$cnames, ctx$rnames), c(ctx$cnames, ctx$rnames))

plot_rel = file_to_tercen(plot_file) %>%
  as_relation() %>%
  as_join_operator(list(), list())

param_rel = file_to_tercen(param_file) %>%
  as_relation() %>%
  as_join_operator(list(), list())

save_relation(list(main_rel, plot_rel, param_rel), ctx)
