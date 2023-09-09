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
#'    the spots are placed using either on their original locations ("spatial")
#'    or user-specified graph layout algorithm.  Default is ("kk" =
#'    [igraph::layout_with_kk()]) algorithm.
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
#' @param edge_color_range,node_color_range bound given in c(min, max) format.
#'   if specified, the ranges are set using these values, instead of the ones
#'   derived from data.  Useful when multiple plots are combined
#'   later for comparison.
#' @param edge_width,node_size specify edge width and node size.  Currently,
#'   default values do not work well.  Might need some manual adjustments.
#' @param edge_alpha,node_alpha specify transparency of edges and nodes
#' @param clip if TRUE (dafault), tissue_img will be cliped to show only
#'   active cell-cell communications
#' @param tissue_img Raster image, like imgRaster(spe)
#' @param image_alpha alpha value for tissue_img. If NA, keep the current value.
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
#' @import ggraph
plot_spatial_ccc_graph <-
  function(ccc_graph,
           graph_layout = "auto",
           cells_of_interest = NULL,
           edges_expanded_to_group = TRUE,
           edge_color = "group_n_edges",
           node_color = "group_n_nodes",
           edge_color_range = NULL,
           node_color_range = NULL,
           edge_width = 0.5,
           edge_alpha = 1,
           node_size = 0.5,
           node_alpha = 1,
           clip = TRUE,
           tissue_img = NULL,
           image_alpha = NA,
           which_on_top = "edge",
           show_arrow = FALSE) {

    if (is.null(ccc_graph)) {
      return(ggplot() + theme_void())
    }


    cells_of_interest_given <- !is.null(cells_of_interest)

    if (cells_of_interest_given) {
      ccc_graph <-
        ccc_graph %>%
        tag_cells_in_ccc_graph(COIs = cells_of_interest,
                               edges_expanded_to_group = edges_expanded_to_group)
    } else {
      ccc_graph <-
        ccc_graph %>%
        tidygraph::activate("nodes") %>%
        tidygraph::mutate(tagged = TRUE) %>%
        tidygraph::activate("edges") %>%
        tidygraph::mutate(tagged = TRUE)
    }

    ccc_graph_nodes <-
      ccc_graph %>%
      tidygraph::activate("nodes") %>%
      tibble::as_tibble()

    ccc_palette_gradient <-
      # RdYlBu is an alternative to Spectral
      grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral")))

    ccc_palette_discrete <- scales::hue_pal()

    # check if node_color is discrete
    node_color_values <-
      ccc_graph %>%
      tidygraph::activate("nodes") %>%
      tidygraph::pull(get(node_color))

    # previously if all node_color_values were integer,
    #   it's considered "discrete";
    # node_color_is_discrete <-
    #   is.factor(node_color_values) |
    #   all(node_color_values == as.integer(node_color_values))

    # if (node_color_is_discrete) {
    #   node_color_is_discrete <-
    #     length(unique(node_color_values)) <= 15
    # }

    # now, simplified; it should be given as factor; otherwise continuous
    node_color_is_factor <- is.factor(node_color_values)
    if (!node_color_is_factor) {
      if (is.null(node_color_range)) {
        max_node_color <- max(node_color_values)
        min_node_color <- min(node_color_values)
      }

      if (max_node_color * min_node_color < 0) {
        mm <- max(abs(max_node_color), abs(min_node_color))
        max_node_color <- mm
        min_node_color <- -mm
      }
    }

    # check if node_color is discrete
    edge_color_values <-
      ccc_graph %>%
      tidygraph::activate("edges") %>%
      tidygraph::pull(get(edge_color))

    # previously if all node_color_values were integer,
    #   it's considered "discrete";
    # node_color_is_discrete <-
    #   is.factor(node_color_values) |
    #   all(node_color_values == as.integer(node_color_values))

    # if (node_color_is_discrete) {
    #   node_color_is_discrete <-
    #     length(unique(node_color_values)) <= 15
    # }

    # now, simplified; it should be given as factor; otherwise continuous
    edge_color_is_factor <- is.factor(edge_color_values)
    if (!edge_color_is_factor) {
      if (is.null(edge_color_range)) {
        max_edge_color <- max(edge_color_values)
        min_edge_color <- min(edge_color_values)
      }

      if (max_edge_color * min_edge_color < 0) {
        mm <- max(abs(max_edge_color), abs(min_edge_color))
        max_edge_color <- mm
        min_edge_color <- -mm
      }
    }
    # node_size_original <- node_size
    # node_alpha_original <- node_alpha
    #
    # edge_width_original <- edge_width
    # edge_alpha_original <- edge_alpha

    edge_alpha_range <- c(0.1, edge_alpha)
    node_alpha_range <- c(0.1, node_alpha)


    if (!is.null(tissue_img)) {
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

      # spatial_image_grob <- grid::rasterGrob(tissue_img)
      spatial_image_grob <-
        grid::rasterGrob(modify_alpha_image(tissue_img,
                                            image.alpha = image_alpha))

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

      # add tissue image
      # when ghost_img was used, this block is used only
      #   when ghost_img == FALSE
      ggraph_ccc <-
        ggraph_ccc +
        geom_spatial(
          data = tibble::tibble_row(sample = "sample",
                                    grob = spatial_image_grob),
          ggplot2::aes(grob = grob),
          x = 0.5,
          y = 0.5
        )
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
          ggraph::ggraph() + scale_y_reverse()
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
      gg <-
        gg +
        ggraph::geom_node_point(
          ggplot2::aes(fill = get(node_color),
              alpha = as.numeric(tagged)),
          shape = 21,
          color = "black",
          size = node_size
        )

      if (node_color_is_factor) {
        unique_node_color_values <- sort(unique(node_color_values))
        node_color_discrete <-
          ccc_palette_discrete(n = length(unique_node_color_values))
        names(node_color_discrete) <- unique_node_color_values

        gg <-
          gg +
          ggplot2::scale_color_manual(name = node_color, values = node_color_discrete) +
          ggplot2::scale_fill_manual(name = node_color, values = node_color_discrete)
      } else {
        gg <-
          gg +
          ggplot2::scale_color_gradientn(
            name = node_color,
            colours = ccc_palette_gradient(n = 15),
            limits = c(min_node_color, max_node_color)
          ) +
          ggplot2::scale_fill_gradientn(
            name = node_color,
            colours = ccc_palette_gradient(n = 15),
            limits = c(min_node_color, max_node_color)
          )
      }

      gg +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = node_size))) +
        ggplot2::scale_alpha(range = node_alpha_range, guide = FALSE) +
        ggplot2::scale_size(guide = FALSE)
    }

    add_ggraph_edges <- function(gg) {
      gg <-
        gg +
        ggraph::geom_edge_link(
          ggplot2::aes(color = get(edge_color),
              alpha = as.numeric(tagged)),

          # for all edges
          # alpha = edge_alpha,
          # width = edge_width,

          # add arrow
          arrow = grid::arrow(
            angle = 15,
            type = "closed",
            length = grid::unit(ifelse(show_arrow,
                                       0.01 * edge_width,
                                       0), "npc")
          )
        )

      if (edge_color_is_factor) {
        unique_edge_color_values <- sort(unique(edge_color_values))
        edge_color_discrete <-
          ccc_palette_discrete(n = length(unique_edge_color_values))
        names(edge_color_discrete) <- unique_edge_color_values

        gg <-
          gg +
          ggraph::scale_edge_color_manual(
            name = edge_color,
            values = edge_color_discrete
          )
      } else {
        gg <-
          gg +
          ggraph::scale_edge_color_gradientn(
            name = edge_color,
            colours = ccc_palette_gradient(n = 15),
            limits = c(min_edge_color, max_edge_color)
          )
      }

      # return ggplot object
      gg + ggraph::scale_edge_alpha(range = edge_alpha_range, guide = FALSE)
    }

    if (which_on_top == "edge") {
      ggraph_ccc %<>%
        add_ggraph_nodes() %>%
        add_ggraph_edges
    } else {
      ggraph_ccc %<>%
        add_ggraph_edges() %>%
        add_ggraph_nodes()
    }

    if (!is.null(tissue_img)) {
      ggraph_ccc <-
        ggraph_ccc +
        # set boundary
        ggplot2::xlim(x_min, x_max) +
        ggplot2::ylim(y_max, y_min) +
        ggplot2::coord_fixed(expand = FALSE) +
        # remove background
        ggplot2::theme_void()
    } else {
      ggraph_ccc <-
        ggraph_ccc +
        # theme_graph()
        # remove background
        ggraph::theme_graph()

      if (graph_layout == "spatial") {
        ggraph_ccc <-
          ggraph_ccc +
          ggplot2::coord_fixed(expand = FALSE)
      }
    }

    ggraph_ccc
  }


#' Plot spatial feature
#'
#' Visualize spatial features
#'   w/ or w/o tissue image
#'
#'   The spots are placed on their original locations,
#'   keeping spatial relations.
#'
#' @param spe SpatialExperiment object
#' @param feature string, column name for a feature to be visualized
#' @param cells_of_interest array of cell ids to highlight.
#' @param spot_size specify size of spot (feature).  Currently,
#'   default values do not work well.  Might need some manual adjustments.
#' @param spot_alpha specify transparency of spot (feature)
#' @param image_alpha alpha value for tissue image.
#'   If NA, keep the current value.
#' @param show_tissue_image if TRUE, tissue image is shown.
#' @param clip if TRUE (dafault), tissue image will be cliped to the size to
#'   cover only the spots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_spatial_features(spe,
#'                       feature = "mclust_6")
#'
#' plot_spatial_features(spe,
#'                       feature = "mclust_6",
#'                       image_alpha = 0.75)
#' }
#'
plot_spatial_feature <-
  function(spe,
           feature,
           cells_of_interest = NULL,
           spot_size = 1.6,
           spot_alpha = 1,
           image_alpha = NA,
           show_tissue_image = TRUE,
           clip = TRUE) {
    spatial_col_data <-
      get_spatial_data(spe)

    cells_of_interest_given <- !is.null(cells_of_interest)

    if (cells_of_interest_given) {
      spatial_col_data <-
        spatial_col_data %>%
        dplyr::filter(cell_id %in% cells_of_interest)
    }

    gp_spatial <-
      spatial_col_data %>%
      ggplot2::ggplot(ggplot2::aes(
        x = spot_x,
        y = spot_y,
        fill = get(feature)
      )) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = feature))

    tissue_img <- SpatialExperiment::imgRaster(spe)

    if (!is.null(tissue_img) & show_tissue_image) {
      img_width <- ncol(tissue_img)
      img_height <- nrow(tissue_img)

      x_min <- 0
      x_max <- img_width
      y_min <- 0
      y_max <- img_height

      if (clip) {
        x_min <- max(x_min, floor(min(spatial_col_data$spot_x)) - 1)
        x_max <-
          min(x_max, ceiling(max(spatial_col_data$spot_x)) + 1)
        y_min <- max(y_min, floor(min(spatial_col_data$spot_y)) - 1)
        y_max <-
          min(y_max, ceiling(max(spatial_col_data$spot_y)) + 1)

        margin <- round(0.05 * min(x_max - x_min, y_max - y_min))
        x_min <- max(x_min - margin, 0)
        x_max <- min(x_max + margin, img_width)
        y_min <- max(y_min - margin, 0)
        y_max <- min(y_max + margin, img_height)
      }

      tissue_img <- tissue_img[y_min:y_max, x_min:x_max]

      # make tissue_img a bit less opaque, mainly aesthetic reason.
      spatial_image_grob <-
        grid::rasterGrob(modify_alpha_image(tissue_img,
                                            image.alpha = image_alpha))
      gp_spatial <-
        gp_spatial +
        geom_spatial(
          data = tibble::tibble_row(sample = "sample",
                                    grob = spatial_image_grob),
          ggplot2::aes(grob = grob),
          x = 0.5,
          y = 0.5
        )
    } else {
      x_min <- floor(min(spatial_col_data$spot_x)) - 1
      x_max <- ceiling(max(spatial_col_data$spot_x)) + 1

      y_min <- floor(min(spatial_col_data$spot_y)) - 1
      y_max <- ceiling(max(spatial_col_data$spot_y)) + 1

      margin <- round(0.05 * min(x_max - x_min, y_max - y_min))
      x_min <- x_min - margin
      x_max <- x_max + margin
      y_min <- y_min - margin
      y_max <- y_max + margin
    }

    gp_spatial +
      ggplot2::geom_point(
        shape = 21,
        color = "grey15",
        size = spot_size,
        alpha = spot_alpha
      ) +

      # set boundary
      ggplot2::xlim(x_min, x_max) +
      ggplot2::ylim(y_max, y_min) +

      ggplot2::coord_fixed(expand = FALSE) +

      # remove background
      ggplot2::theme_void() +

      scale_fill_continuous(type = "viridis") +
      scale_color_continuous(type = "viridis")
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
#' @importFrom tibble tibble
#'
#' @noRd
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

  GeomCustom <- ggplot2::ggproto(
    "GeomCustom",
    ggplot2::Geom,
    setup_data = function(self, data, params) {
      data <- ggplot2::ggproto_parent(ggplot2::Geom, self)$setup_data(data, params)
      data
    },

    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x = data$x, y = data$y)
      g <- grid::editGrob(data$grob[[1]], vp = vp)
      ggname("geom_spatial", g)
    },

    required_aes = c("grob", "x", "y")
  )

  ggplot2::layer(
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

#' Modify image with alpha
#'
#' @param raster Raster image, like imgRaster(spe)
#' @param image.alpha alpha value for tissue_img. If NA, keep the current value.
#'
#' @noRd
modify_alpha_image <- function(raster, image.alpha = NA) {
  matrix(data = scales::alpha(raster, alpha = image.alpha),
         nrow = nrow(x = raster),
         ncol = ncol(x = raster),
         byrow = TRUE)
}



