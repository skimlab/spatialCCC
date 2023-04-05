#' Plot spatial CCC graph
#'
#' Visualize spatial CCC graph in various formats,
#'   w/ or w/o tissue image
#'
#'   * CCC graph with tissue image. In this format,
#'   the spots are placed on their original locations, keeping
#'   spatial relations.
#'
#'   * CCC graph without tissue image. In this format,
#'    the spots are placed using user-specified graph layout
#'    algorithm.  Default is ("kk")[igraph::layout_with_kk()] algorithm.
#'    See [create_layout()] for available algorithms.
#'
#' @param ccc_graph spatial CCC graph for a LR pair
#' @param graph_layout string for graph layout algorithm. Default is "auto";
#'   see Details for more information
#' @param cells_of_interest array of cell ids to highlight.
#' @param edges_expanded_to_group if TRUE and cells_of_interest is given,
#'   the nodes and the edges are expanded to include all the nodes and edges that
#'   belong to the groups that cells_of_interest belong to.
#' @param edge_color,node_color column name to be used
#'   for edge and node colors, respectively
#' @param edge_range,node_range bound given in c(min, max) format.
#'   if specified, the ranges are set using these values, instead of the ones
#'   derived from data.  Useful when multiple plots are combined
#'   later for comparison.
#' @param edge_width,node_size specify edge width and node size.  Currently,
#'   default values do not work well.  Might need some manual adjustments.
#' @param edge_alpha,node_alpha specify transparency of edges and nodes
#' @param clip if TRUE, tissue_img will be cliped to show only active cell-cell
#'   communications
#' @param tissue_img Raster image, like imgRaster(spe)
#' @param ghost_img if TRUE, tissue_img is not displayed, only used for
#'   placing spots.
#' @param which_on_top either "edge" or "node"
#' @param show_arrow if TRUE, arrow is added to edge.  Default is FALSE.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[["App_Dcc"]],
#'                        tissue_img = imgRaster(spe),
#'                        node_color = "group_diameter",
#'                        node_size = 1,
#'                        edge_color = "group_diameter",
#'                        which_on_top = "edge")
#'
#' plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[["App_Dcc"]],
#'                        node_color = "group_diameter",
#'                        node_size = 1,
#'                        edge_color = "group_diameter",
#'                        which_on_top = "edge")
#'
#' plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[["App_Dcc"]],
#'                        graph_layout = "stress",
#'                        node_color = "group_diameter",
#'                        node_size = 1,
#'                        edge_color = "group_diameter",
#'                        which_on_top = "edge")
#' }
#'
plot_spatial_ccc_graph <-
  function(ccc_graph,
           graph_layout = "auto",
           cells_of_interest = NULL,
           edges_expanded_to_group = TRUE,
           edge_color = "group_n_edges",
           node_color = "group_n_nodes",
           edge_range = NULL,
           node_range = NULL,
           edge_width = 0.5,
           edge_alpha = 1,
           node_size = 0.5,
           node_alpha = 1,
           clip = FALSE,
           tissue_img = NULL,
           ghost_img = FALSE,
           which_on_top = "edge",
           show_arrow = FALSE) {
    cells_of_interest_given <- !is.null(cells_of_interest)

    if (cells_of_interest_given) {
      ccc_graph %<>%
        tidygraph::activate(nodes) %>%
        tidygraph::mutate(InCOI = name %in% cells_of_interest) %>%
        tidygraph::mutate(InFocus = InCOI) %>%
        tidygraph::activate(edges) %>%
        tidygraph::mutate(InCOI = src %in% cells_of_interest |
                            dst %in% cells_of_interest) %>%
        tidygraph::mutate(InFocus = InCOI)

      if (edges_expanded_to_group) {
        groups_of_interest <-
          ccc_graph %>%
          tidygraph::activate(nodes) %>%
          tidygraph::filter(InFocus) %>%
          tidygraph::pull(group) %>%
          unique()

        cells_in_GOI <-
          ccc_graph %>%
          tidygraph::activate(nodes) %>%
          tidygraph::filter(group %in% groups_of_interest) %>%
          tidygraph::pull(name) %>%
          unique()

        ccc_graph %<>%
          tidygraph::activate(edges) %>%
          tidygraph::mutate(InGOI = src %in% cells_in_GOI |
                              dst %in% cells_in_GOI) %>%
          tidygraph::mutate(InFocus = InGOI)
      }
    } else {
      ccc_graph %<>%
        tidygraph::activate(nodes) %>%
        tidygraph::mutate(InFocus = TRUE) %>%
        tidygraph::activate(edges) %>%
        tidygraph::mutate(InFocus = TRUE)
    }

    # ccc_graph %<>%
    #   activate(nodes) %>%
    #   filter(InFocus)

    ccc_graph_nodes <-
      ccc_graph %>%
      tidygraph::activate(nodes) %>%
      tibble::as_tibble()

    ccc_edge_palette <-
      # RdYlBu is an alternative to Spectral
      grDevices::colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))

    if (is.null(edge_range)) {
      edge_color_values <-
        ccc_graph %>%
        tidygraph::activate(edges) %>%
        tidygraph::pull(get(edge_color))

      min_edge_color <- min(edge_color_values)
      max_edge_color <- max(edge_color_values)
    }

    if (max_edge_color * min_edge_color < 0) {
      mm <- max(abs(max_edge_color), abs(min_edge_color))
      max_edge_color <- mm
      min_edge_color <- -mm
    }

    ccc_node_palette <-
      grDevices::colorRampPalette(rev(brewer.pal(7, "Spectral")))

    # check if node_color is discrete
    node_color_values <-
      ccc_graph %>%
      tidygraph::activate(nodes) %>%
      tidygraph::pull(get(node_color))

    node_color_is_discrete <-
      all(node_color_values == as.integer(node_color_values))

    if (node_color_is_discrete) {
      node_color_is_discrete <-
        length(unique(node_color_values)) <= 15
    }

    if (is.null(node_range)) {
      max_node_color <- max(node_color_values)
      min_node_color <- min(node_color_values)
    }

    if (max_node_color * min_node_color < 0) {
      mm <- max(abs(max_node_color), abs(min_node_color))
      max_node_color <- mm
      min_node_color <- -mm
    }


    node_size_original <- node_size
    node_alpha_original <- node_alpha

    edge_width_original <- edge_width
    edge_alpha_original <- edge_alpha

    if (which_on_top == "edge") {
      node_alpha = 0.5 * node_alpha_original
      #node_size = 0.5 * node_size
    } else {
      edge_alpha = 0.5 * edge_alpha_original
      #edge_width = 0.5 * edge_width
    }

    edge_alpha_range <- c(0.1, edge_alpha)
    node_alpha_range <- c(0.1, node_alpha)


    if (!is.null(tissue_img)) {
      #
      # add tissue image
      #
      # tissue_img <- imgRaster(spe)
      img_width <- ncol(tissue_img)
      img_height <- nrow(tissue_img)

      x_min <- 0
      x_max <- img_width
      y_min <- 0
      y_max <- img_height

      if (clip) {
        x_min <- max(x_min, floor(min(ccc_graph_nodes$spot_x)) - 1)
        x_max <-
          min(x_max, ceiling(max(ccc_graph_nodes$spot_x)) + 1)
        y_min <- max(y_min, floor(min(ccc_graph_nodes$spot_y)) - 1)
        y_max <-
          min(y_max, ceiling(max(ccc_graph_nodes$spot_y)) + 1)

        margin <- round(0.05 * min(x_max - x_min, y_max - y_min))
        x_min <- max(x_min - margin, 0)
        x_max <- min(x_max + margin, img_width)
        y_min <- max(y_min - margin, 0)
        y_max <- min(y_max + margin, img_height)

      }

      tissue_img <- tissue_img[y_min:y_max, x_min:x_max]

      spatial_image_grob <- grid::rasterGrob(tissue_img)

      #
      # create ggraph object
      #
      graph_layout <- "manual"
      ggraph_ccc <-
        ggraph::create_layout(
          graph = ccc_graph,
          layout = "manual",
          ## ignore graph_layout option
          x = spot_x,
          y = spot_y
        ) %>%
        ggraph::ggraph()

      if (!ghost_img) {
        # add tissue image
        ggraph_ccc <-
          ggraph_ccc +
          geom_spatial(
            data = tibble::tibble_row(sample = "sample",
                                      grob = spatial_image_grob),
            aes(grob = grob),
            x = 0.5,
            y = 0.5
          )
      }
    } else {
      #
      # create ggraph object
      #
      if (graph_layout == "spatial") {
        ggraph_ccc <-
          ggraph::create_layout(
            graph = ccc_graph,
            layout = "manual",
            ## ignore graph_layout option
            x = spot_x,
            y = spot_y
          ) %>%
          ggraph::ggraph()
      } else {
        # some layout alg. does not work
        # so, default back to stress
        # otherwise, left to the user
        if (graph_layout == "auto")
          graph_layout <- "kk"

        ggraph_ccc <-
          ggraph::ggraph(graph = ccc_graph,
                         layout = graph_layout)
      }
    }

    add_ggraph_nodes <- function(gg) {
      if (node_color_is_discrete) {
        unique_node_color_values <- sort(unique(node_color_values))
        node_color_discrete <-
          ccc_node_palette(n = length(unique_node_color_values))
        names(node_color_discrete) <- unique_node_color_values

        gg +
          ggraph::geom_node_point(
            aes(# color = factor(get(node_color)),
              fill = factor(get(node_color)),
              alpha = as.numeric(InFocus)),
            pch = 21,
            color = "black",
            size = node_size
          ) +
          scale_color_manual(name = "node", values = node_color_discrete) +
          scale_fill_manual(name = "node", values = node_color_discrete) +
          guides(color = guide_legend(override.aes = list(size = node_size_original * 2))) +
          scale_alpha(guide = FALSE) +
          scale_size(guide = FALSE)
      } else {
        gg +
          ggraph::geom_node_point(
            aes(
              # color = get(node_color),
              fill = get(node_color),
              alpha = as.numeric(InFocus)
            ),
            pch = 21,
            color = "black",
            size = node_size
          ) +

          # for all nodes
          # alpha = node_alpha,
          # size = node_size
          scale_color_gradientn(
            name = "node",
            colours = ccc_node_palette(n = 15),
            limits = c(min_node_color, max_node_color)
          ) +
          scale_fill_gradientn(
            name = "node",
            colours = ccc_node_palette(n = 15),
            limits = c(min_node_color, max_node_color)
          ) +
          guides(color = guide_legend(override.aes = list(size = node_size_original * 2))) +
          scale_alpha(range = node_alpha_range, guide = FALSE) +
          scale_size(guide = FALSE)
      }
    }

    add_ggraph_edges <- function(gg) {
      gg +
        ggraph::geom_edge_link(
          aes(color = get(edge_color),
              alpha = as.numeric(InFocus)),

          # for all edges
          # alpha = edge_alpha,
          # width = edge_width,

          # add arrow
          arrow = grid::arrow(
            angle = 15,
            type = "closed",
            length = unit(ifelse(show_arrow,
                                 0.01 * edge_width,
                                 0), "npc")
          )
        ) +
        ggraph::scale_edge_color_gradientn(
          name = "edge",
          colours = ccc_edge_palette(n = 15),
          limits = c(min_edge_color, max_edge_color)
        ) +
        ggraph::scale_edge_alpha(range = edge_alpha_range, guide = FALSE)
    }

    if (which_on_top == "edge") {
      ggraph_ccc %<>%
        add_ggraph_nodes() %>%
        add_ggraph_edges
      # guides(size = FALSE, alpha = FALSE)
    } else {
      ggraph_ccc %<>%
        add_ggraph_edges() %>%
        add_ggraph_nodes()
    }

    if (!is.null(tissue_img)) {
      ggraph_ccc <-
        ggraph_ccc +
        # set boundary
        xlim(x_min, x_max) +
        ylim(y_max, y_min) +
        coord_fixed(expand = FALSE) +
        # remove background
        ggplot2::theme_void()
    } else {
      ggraph_ccc <-
        ggraph_ccc +
        # theme_graph()
        # remove background
        theme_graph()
    }

    ggraph_ccc
  }

#' A ggplot2 layer for visualizing the Visium histology
#'
#' This function defines a [ggplot2::layer()] for visualizing the histology
#' image from Visium. It can be combined with other ggplot2 functions for
#' visualizing the clusters as in [plot_spatial_ccc_graph()].
#'
#' @param mapping Passed to `ggplot2::layer(mapping)` where `grob`, `x` and `y`
#' are required.
#' @param data Passed to `ggplot2::layer(data)`.
#' @param stat Passed to `ggplot2::layer(stat)`.
#' @param position Passed to `ggplot2::layer(position)`.
#' @param na.rm Passed to `ggplot2::layer(params = list(na.rm))`.
#' @param show.legend Passed to `ggplot2::layer(show.legend)`.
#' @param inherit.aes Passed to `ggplot2::layer(inherit.aes)`.
#' @param ... Other arguments passed to `ggplot2::layer(params = list(...))`.
#'
#' @return A [ggplot2::layer()] for the histology information.
#' @author 10x Genomics
#' @export
#' @importFrom tibble tibble
#'
geom_spatial <- function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  ## To avoid a NOTE on R CMD check
  ggname <- function(prefix, grob) {
    grob$name <- grid::grobName(grob, prefix)
    grob
  }

  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },

    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x = data$x, y = data$y)
      g <- grid::editGrob(data$grob[[1]], vp = vp)
      ggname("geom_spatial", g)
    },

    required_aes = c("grob", "x", "y")
  )

  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
