#' @import Seurat
#' @import tidyverse
#' @import networkD3

#' Create a Sankey plot
#' 
#' @param integ_srat Seurat dataset objects
#' @param clusters String denoting column to get Sankey plot of
#' @param out_directory String of output directory


sankey_plot <- function(
  integ_srat,
  clusters = "seurat_clusters"
  #out_directory
) {
    integ_human <- subset(integ_srat, species == "human")
    integ_mouse <- subset(integ_srat, species == "mouse")
    
    unique_cell_types_human <- unique(na.omit(c(integ_human$common_cell_name)))
    unique_cell_types_mouse <- unique(na.omit(c(integ_mouse$common_cell_name)))
    unique_cell_types_mouse <- paste("Mouse:", unique_cell_types_mouse, sep = " ")

    # determining name of unique clusters
    unique_cell_clusters <- unique(as.character(c(integ_srat@meta.data[[clusters]])))
    unique_cell_clusters <- paste("Cluster", unique_cell_clusters, sep = " ")

    # appending all names into a data frame
    nodes <- data.frame("name" = 
        c(unique_cell_types_human, # All human cell types
        unique_cell_types_mouse, # All integrated clusters
        unique_cell_clusters))# All mouse cell types

    # creating links between human cells to cluster
    human_to_cluster_links <- as.data.frame(cbind(
        as.character(integ_human$common_cell_name),
        paste0("Cluster ", as.character(integ_human@meta.data[[clusters]]))
    ))

    # creating links between cluster and mice cells
    cluster_to_mouse_links <- as.data.frame(cbind(
        paste0("Cluster ", as.character(integ_mouse@meta.data[[clusters]])),
        paste0("Mouse: ",as.character(integ_mouse$common_cell_name)) 
    ))

    # appending all links and removing links with NA
    links <- rbind(human_to_cluster_links, cluster_to_mouse_links, na.rm = TRUE) %>%
        table() %>%
        as.data.frame() %>%
        # filter paths based on minimum cell count
        subset(Freq > 200) %>%
        lapply(as.character)

    names(links) = c("source", "target", "value")

    # determining which nodes are left after filtering 
    all_nodes <- unique(c(links$source, links$target))
    used_nodes <- unique(c(links$source[links$value > 0], links$target[links$value > 0]))
    used_nodes <- data.frame("name" = used_nodes)

    # converting the nodes into their index in the data frame
    # Note: sankeyNetwork requires 0 indexing, and data frames index starting at 1, therefore needed to minus 1 to the value
    source_num <- match(links$source, used_nodes$name) - 1
    target_num <- match(links$target, used_nodes$name) - 1

    # updating links from name to number
    links$source <- source_num
    links$target <- target_num

    links <- as.data.frame(links)

    # creating sankey plot
    sankey <- sankeyNetwork(Links = links, Nodes = used_nodes,
    Source = "source", Target = "target",
    Value = "value", NodeID = "name",
    fontSize= 10, nodeWidth = 20,
    nodePadding = 10)

    # generates the number of cells per path that appear on the plot
    library(htmlwidgets)
    htmlwidgets::onRender(sankey, '
        function(el) { 
        var nodeWidth = this.sankey.nodeWidth();
        var links = this.sankey.links();
        var nodes = this.sankey.nodes();

        links.forEach((d, i) => {
            if (d.value > 400) {
            var startX = d.source.x + nodeWidth;
            var endX = d.target.x;
            
            var startY = d.source.y + d.sy + d.dy / 2;
            var endY = d.target.y + d.ty + d.dy / 2;
            
            d3.select(el).select("svg g")
                .append("text")
                .attr("text-anchor", "middle")
                .attr("alignment-baseline", "middle")
                .attr("x", startX + ((endX - startX) / 2))
                .attr("y", startY + ((endY - startY) / 2))
                .attr("font-size", "10px")
                .text(d.value);
            }
        });
        }
    ')
    #html_file <- file.path(out_directory, "sankey_temp.html")
    #saveWidget(sankey, html_file)
    # webshot(html_file, file = paste(out_directory, "/sankey_diagram.pdf", sep = ""))
}


