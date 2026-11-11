# Calibration & Quantification app

Web app for building calibration curves (linear, quadratic, log, power; OLS or
1/X, 1/X² weighting), quantifying samples and comparing groups.

**Use it online:** https://YOUR-GITHUB-USERNAME.github.io/calibration-app/

- No installation needed: R runs inside your web browser (Shinylive / webR).
- The first load takes 20-60 s (R is downloaded into the browser), later visits are faster.
- **Your data never leaves your computer**: files are processed locally in the browser.
- Click *Load example data* in tab 1 to try it.

## Input files

CSV files with `;`, `,` or tab separators; decimal commas are accepted (French Excel exports).

| File | Required columns | Optional columns |
|------|------------------|------------------|
| Calibration | concentration, response (area/signal) | - |
| Samples | sample name, response | `Theoretical_Concentration`, `Dilution_Factor`, any grouping column (treatment, class...) |

Templates can be downloaded from tab 2.

## For maintainers

- `app/app.R` is the application (plus example data in `app/`).
- `docs/` is the generated website served by GitHub Pages; do not edit it by hand.
- After editing `app/app.R`: `Rscript build.R`, check locally with
  `Rscript -e "httpuv::runStaticServer('docs')"`, then commit and push.
- The app also runs in desktop R: `shiny::runApp("app")`.
