# Pure helper functions of the ComBat operator (no Tercen calls), so they can be
# used by main.R and by the dev/ scripts alike.

# Split "a; b" into c("a", "b"); spaces around the separator are ignored.
split_values = function(x, sep = ";") {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, sep, fixed = TRUE)[[1]])
}

# "a, b, c" for messages; long lists are cut after `max` items.
list_text = function(x, max = 15) {
  if (length(x) <= max) return(paste(x, collapse = ", "))
  paste0(paste(x[1:max], collapse = ", "), ", ... (", length(x) - max, " more)")
}

# One batch label per sample: the colour factor values joined with ";".
batch_labels = function(colour_df) {
  do.call(paste, c(lapply(colour_df, function(v) trimws(as.character(v))), sep = ";"))
}

# Peptide x sample matrix from the long crosstab table (.ri, .ci zero-based).
# Row and column order follow .ri and .ci, as acast() did in the Shiny operator.
build_matrix = function(ri, ci, y, n_rows, n_cols) {
  if (anyDuplicated(cbind(ri, ci))) {
    stop("More than one value per peptide and sample. Check that the column factors ",
         "identify one sample per column (e.g. Barcode and Array/Row).", call. = FALSE)
  }
  X = matrix(NA_real_, n_rows, n_cols, dimnames = list(rowSeq = 0:(n_rows - 1), colSeq = 0:(n_cols - 1)))
  X[cbind(ri + 1, ci + 1)] = y
  X
}

# The batch of each sample (column). All cells of a column must share one batch.
sample_batches = function(ci, labels, n_cols) {
  per_col = tapply(labels, factor(ci, levels = 0:(n_cols - 1)), unique, simplify = FALSE)
  multi = which(lengths(per_col) > 1)
  if (length(multi) > 0) {
    stop("A sample (crosstab column) belongs to more than one batch. The batch factors on ",
         "colour must have one value per column; put the factors that define the sample on the columns.",
         call. = FALSE)
  }
  vapply(per_col, function(v) if (length(v)) v else NA_character_, "")
}

# Resolve ReferenceBatch against the batch labels; stops with a hint to valid batches.
resolve_reference = function(ref_string, batches, colour_names) {
  given = split_values(ref_string)
  want = paste(given, collapse = ";")
  available = sort(unique(batches))
  if (nzchar(want) && want %in% available) return(want)
  n_col = length(colour_names)
  if (n_col > 1) {
    how = paste0("With several batch factors (", paste(colour_names, collapse = ", "),
                 ") give one value per factor in this order, separated by ';'. ",
                 "An example available batch: ", available[1], ".")
    if (length(given) != n_col) {
      stop(if (length(given) == 1) "Only 1 value was" else paste(length(given), "values were"),
           " given as reference batch, but ", n_col, " colour factors were used. ", how, call. = FALSE)
    }
    stop("ReferenceBatch '", ref_string, "' is not one of the batches. ", how, call. = FALSE)
  }
  stop("ReferenceBatch '", ref_string, "' is not one of the batches. Available batches: ",
       list_text(available), call. = FALSE)
}

# Find the column factor for a user-given name, with or without namespace prefix.
match_factor = function(name, cnames) {
  hit = cnames[cnames == name]
  if (length(hit) == 0) hit = cnames[endsWith(cnames, paste0(".", name))]
  if (length(hit) == 1) return(hit)
  if (length(hit) == 0) {
    stop("The fit condition uses the factor '", name, "', but it is not on the crosstab columns. ",
         "Add '", name, "' to the columns, or correct the FitConditionFactors setting.", call. = FALSE)
  }
  stop("The fit-condition factor '", name, "' matches several column factors (",
       paste(hit, collapse = ", "), "). Use the full name in FitConditionFactors.", call. = FALSE)
}

# Logical vector over samples: TRUE for samples of the fit condition.
fit_condition_samples = function(factor_string, value_string, col_df) {
  names_wanted = split_values(factor_string)
  if (length(names_wanted) == 0) {
    stop("UseFitCondition is on, but FitConditionFactors is empty. Give the column factor(s) ",
         "that define the fit condition, e.g. Grouping.", call. = FALSE)
  }
  cols = vapply(names_wanted, match_factor, "", cnames = names(col_df))
  combos = split_values(value_string, "|")
  if (length(combos) == 0) {
    stop("UseFitCondition is on, but FitConditionValues is empty. Give the value(s) of ",
         paste(names_wanted, collapse = ";"), " for the fit condition, e.g. REF.", call. = FALSE)
  }
  sample_combo = batch_labels(col_df[, cols, drop = FALSE])
  wanted = vapply(combos, function(cmb) {
    v = split_values(cmb)
    if (length(v) != length(cols)) {
      stop("The fit condition '", cmb, "' has ", length(v), " value(s), but FitConditionFactors names ",
           length(cols), " factor(s) (", paste(names_wanted, collapse = ";"), "). Give one value per ",
           "factor, separated by ';'.", call. = FALSE)
    }
    paste(v, collapse = ";")
  }, "")
  sel = sample_combo %in% wanted
  if (!any(sel)) {
    stop("No sample matches the fit condition ", paste(wanted, collapse = " | "), ". Values of ",
         paste(names_wanted, collapse = ";"), " in the data: ", list_text(sort(unique(sample_combo))),
         call. = FALSE)
  }
  sel
}

# Plain-language check that every batch has enough fit samples for the model.
check_fit_counts = function(fit_batches, all_batches, mean_only, fit_condition) {
  counts = table(factor(fit_batches, levels = sort(unique(all_batches))))
  shown = if (any(counts < 2)) counts[order(counts)] else counts
  listing = list_text(sprintf("%s: %d", names(shown), as.integer(shown)))
  what = if (fit_condition) "fit-condition samples" else "samples"
  if (any(counts == 0)) {
    stop("Some batches have no ", what, ", so they cannot be corrected (", listing, "). ",
         "Every batch needs at least one.", call. = FALSE)
  }
  if (!mean_only && any(counts < 2)) {
    stop("The L/S model needs at least 2 ", what, " in every batch (", listing, "). ",
         "Use the L model, or ", if (fit_condition) "add more conditions to the fit condition." else "remove batches with one sample.",
         call. = FALSE)
  }
  if (mean_only && all(counts < 2)) {
    stop("At least one batch needs 2 or more ", what, " to estimate the normal spread of the data (",
         listing, "). ", if (fit_condition) "Add more conditions to the fit condition." else "",
         call. = FALSE)
  }
  invisible(counts)
}

# The settings as shown to the user, in the PCA plot caption and the parameters CSV.
# Several fit conditions ("a | b") are written as "a or b", so "|" only separates
# the parameters in the caption.
parameter_table = function(mean_only, batch_factors, ref, fit_factors, fit_values) {
  fit = if (is.null(fit_factors)) "not used" else
    paste(paste(split_values(fit_factors), collapse = ";"), "-",
          paste(vapply(split_values(fit_values, "|"), function(v) paste(split_values(v), collapse = ";"), ""),
                collapse = " or "))
  data.frame(parameter = c("Model type", "Batch factors", "Reference batch", "Fit condition"),
             value = c(if (mean_only) "L" else "L/S", paste(batch_factors, collapse = ";"),
                       if (is.null(ref)) "not used" else ref, fit))
}

# One line for the plot caption (without the batch factors, which the colour legend shows).
parameter_text = function(params) {
  shown = params[params$parameter != "Batch factors", ]
  paste("Parameters:", paste(shown$parameter, shown$value, sep = ": ", collapse = " | "))
}

# PC1 vs PC2 before and after correction, coloured by batch, with the settings
# as caption under the figures.
pca_plot = function(X0, Xc, batches, fit_sample = NULL, caption = NULL) {
  pcs = function(X, stage) {
    p = prcomp(t(X))
    data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], batch = batches, stage = stage)
  }
  df = rbind(pcs(X0, "before"), pcs(Xc, "after"))
  df$stage = factor(df$stage, levels = c("before", "after"))
  if (is.null(fit_sample)) {
    plt = ggplot(df, aes(PC1, PC2, colour = batch)) + geom_point(size = 2.5)
  } else {
    df$fit = factor(ifelse(rep(fit_sample, 2), "fit condition", "other"), levels = c("fit condition", "other"))
    plt = ggplot(df, aes(PC1, PC2, colour = batch, shape = fit)) + geom_point(size = 2.5) +
      scale_shape_manual(values = c(`fit condition` = 17, other = 16), name = NULL)
  }
  plt + facet_wrap(~stage, scales = "free") + theme_bw() +
    labs(caption = caption) +
    theme(plot.caption = element_text(hjust = 0, size = 10, margin = margin(t = 10)),
          plot.caption.position = "plot")
}
