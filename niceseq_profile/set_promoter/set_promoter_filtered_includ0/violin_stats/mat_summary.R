rm(list = ls())


library(Rseb)

devtools::install_github("sebastian-gregoricchio/Rseb",
			 build_manual = TRUE,
                         build_vignettes = TRUE)



                        density.profile.by.group <-
  Rseb::plot.density.profile(
    matrix.file = deeptools.matrix,
    signal.type = "mean",
    error.type = "sem",
    plot.by.group = T,
    missing.data.as.zero = T,
    y.identical.auto = T,
    text.size = 10,
    axis.line.width = 0.25,
    line.width = 0.5,
    plot.error = T,
    write.reference.points = T,
    plot.vertical.lines = F,
    colors = c("steelblue", "mediumseagreen"),
    n.row.multiplot = 2,
    by.row = T,
    y.lab = "Mean ChIP-seq signal \u00b1 SEM")