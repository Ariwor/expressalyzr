#' ---
#' title: "Gating report"
#' output: html_document
#' ---

#+ echo=F, warnings=F, message=T, fig.width = 15, fig.height = 15
cat(paste0("#' ", experiment_name))

pops <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")

for (pop in pops[2:length(pops)]) {
  p <- ggcyto::autoplot(gs, pop)
  print(p)
}
