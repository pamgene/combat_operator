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

Pinned smoke test, run by the release install check and the library gate (not the scientific test — that is `dev/compare_old_new.R` and the preview).

- Input `tests/test_in.csv` = `dev/data/example_data_QC.csv` reduced to `Barcode`, `Row`, `ID`, `value`. Projection: rows `ID`, columns `Barcode` + `Row`, colour `Barcode`, y `value`.
- Every setting pinned at its default (`L/S`, no reference batch, no fit condition). `equalityMethod: R2` (0.99).
- Four expected outputs, in this order (the CmbCor table is joined to the crosstab's column and row tables, so these count as outputs): `test_out_1` CmbCor (`.rids`, `.cids`, `ds0.CmbCor`), `test_out_2` column table, `test_out_3` row table, `test_out_4` PCA plot. `filename` of the plot is in `skipColumns` (R temp file name); the PNG bytes are compared (two runs gave identical bytes).
- Each expected CSV has a `.csv.schema` sidecar (column types); doubles are written type-faithfully (`641031403.0`).

### Regenerating the expected output (intended behaviour change, e.g. phase B bug fixes)

1. Push the change; CI builds `ghcr.io/pamgene/combat_operator:<short sha>` for any branch.
2. In Tercen Studio (`tercen-studio` skill): install the operator from that commit without test (`install_operator`, `testRequired: false`), import `tests/test_in.csv`, build a data step with the projection and pinned settings above, run it **twice**.
3. Export the four output relations of the step's computed relation (CSV + schema JSON with `id`/`rev` removed); both runs must be identical.
4. Replace `tests/test_out_*`, commit, then run the real gate: `tercenctl operator install -r https://github.com/pamgene/combat_operator -t <sha> --rm` must end with `default_params_LS successful`.

Studio note (Oct 2026): on the Docker Desktop WSL2 kernel 6.6.87.2 the Studio worker cannot start operator containers (`netavark: nftables error`). Workaround: inside the worker container, `/etc/containers/containers.conf.d/90-host-network.conf` with `[containers]` / `netns="host"`, then `docker restart tercen_studio-tercen-worker-1` (not recreate — that removes the file).

## Release rules

1. Never point `operator.json`'s `container` at `:main`/`:latest`.
2. Before tagging `X.Y.Z`: set `"container": "ghcr.io/pamgene/combat_operator:X.Y.Z"`, commit, then tag `X.Y.Z`.
3. The release install-check is fatal by design; a red release burns a patch number.
4. Push to `main`/`master`: CI builds the image tagged with branch name and full/short commit SHA (GHA layer cache; dependencies are a separate layer, rebuilt only when `renv.lock` changes).

## Memory model (`memory_model.json`)

First-principles cut: 300 MB baseline + 64 B per crosstab cell (ComBat keeps ~8 copies of the data matrix). Refit from task metas if exit-137 errors appear.
