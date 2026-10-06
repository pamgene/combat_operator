# combat_operator — maintenance guide

Tercen R operator (`main.R`) that removes technical batch effects with ComBat (`R/pgcombat.R`, L or L/S model, optional reference batch, optional fit condition) and returns `CmbCor` plus a before/after PCA PNG. Runs on `tercen/runtime-r44` (R 4.4.3); packages pinned in `renv.lock`. Terms: see `CONTEXT.md`.

## Layout

- `main.R`: Tercen I/O and orchestration only.
- `R/functions.R`: pure helpers (matrix building, settings validation, error texts, PCA plot).
- `R/pgcombat.R`: the ComBat implementation (R6 `pgCombat`, `fit()` / `apply()`).
- `dev/`: test data (`dev/data/`), mock Tercen context, preview and comparison scripts. Not in the image.
- `dev/reference/`: frozen copy of the old Shiny operator code (commit `ecfbdb7`), used by `dev/compare_old_new.R`. Do not edit.

## Checks after every change

```bash
docker build -t combat_operator:dev .
docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/preview.R
docker run --rm -v "$PWD:/src" --entrypoint Rscript combat_operator:dev /src/dev/compare_old_new.R
```

- `preview.R` runs the unmodified `main.R` on all data sets and settings (incl. error cases): look at `dev/preview/index.html` and `dev/preview/errors.txt` and show the PNGs to the user.
- `compare_old_new.R` compares CmbCor with the old Shiny operator. L-model cases must stay equal (≤ 1e-10).
- `dev/data/example_data_conditionfit_ref.csv`: `Sample.name` is recoded (2-letter codes). Never add files with real sample/patient names.

## The Tercen unit test (`tests/test.json`)

Not created yet. When added: pinned smoke test with golden output from a real run (see the create-operator skill).

## Release rules

1. Never point `operator.json`'s `container` at `:main`/`:latest`.
2. Before tagging `X.Y.Z`: set `"container": "ghcr.io/pamgene/combat_operator:X.Y.Z"`, commit, then tag `X.Y.Z`.
3. The release install-check is fatal by design; a red release burns a patch number.
4. Push to `main`/`master`: CI builds the image tagged with branch name and full/short commit SHA (GHA layer cache; dependencies are a separate layer, rebuilt only when `renv.lock` changes).

## Memory model (`memory_model.json`)

First-principles cut: 300 MB baseline + 64 B per crosstab cell (ComBat keeps ~8 copies of the data matrix). Refit from task metas if exit-137 errors appear.
