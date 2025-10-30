plot_elbow.Seurat.elbow.plot <- function(obj, ndim=30,outputPrefix,W=5,H=3) { 

      
	p <- Seurat::ElbowPlot(obj,ndim=ndim);
	p <- p + theme(panel.grid = element_line(color = "gray", linewidth= 0.25, linetype = 1));

	FF <- paste0(outputPrefix,".pca_elbow_plot.pdf"); cat("\nprinting ", FF, "\n");
	pdf(FF, width=W, height=H);print(p); dev.off();

}