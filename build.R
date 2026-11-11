# Rebuild the website after editing app/app.R
# Run from this folder:  Rscript build.R
# Then test locally:     Rscript -e "httpuv::runStaticServer('docs')"
# and publish:           git add -A && git commit -m "Update app" && git push

unlink("docs", recursive = TRUE)
shinylive::export("app", "docs")
