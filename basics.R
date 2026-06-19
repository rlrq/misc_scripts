## this function moves the foremost geom layer to the back
## (i.e. the layer that's rendered last is now rendered first.)
## this is particularly useful when annotation portions of a ggtree tree with rectangles/shapes, as annotation layers are plotted ON TOP of the tree by default and this leads to the tree becoming obscured. Even though the annotation shapes are non-opaque, they still change the colour of the branches that are below them, so I prefer the annotation layers to be under the tree itself.
send_layer_to_back <- function(ggplot_obj){
    num_layers <- length(ggplot_obj$layers)
    reordered_layers <- c(num_layers, 1:(num_layers-1))
    ggplot_obj$layers <- ggplot_obj$layers[reordered_layers]
    return(ggplot_obj)
}

## make path
mkpath <- function(...){paste(..., sep = '/')}

## get file extension
file_ext <- function(fname){
    ext <- stringr::str_extract(fname, "(?<=/.)[^.]+$")
    if (is.na(ext)){return('')}
    return(ext)
}

## function to save plots with default PDF+PNG being default formats
save_plot_noDefaultDir <- function(fname, p, w = 9.91, h = 5.44, units = "in", fmt = NA,
                                   dir = NA, ...){
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    if(is.na(fmt)){
        if (file_ext(fname) != ''){
            print(h)
            ggplot2::ggsave(paste0(dir, "/", fname), 
                            p, w = w, h = h, units = units, ...) ##9.91x5.44in if legend.pos at side
            return()
        } else {
            fmt = c("pdf", "png")
        }
    }
    for (ext in fmt){
        ggplot2::ggsave(paste0(dir, "/", fname, ".", ext), 
                        p, w = w, h = h, units = units, ...) ##9.91x5.44in if legend.pos at side
    }
}

## function to label grobs in top left corner to generate labelled subplots
label_subplot_grob <- function(lab, p, fontsize = 10, x = unit(0.5, "npc"), y = unit(1, "npc")){
    gridExtra::arrangeGrob(grobs = list(plot = p),
                           left = grid::textGrob(rlang::expr(bold(!!lab)),
                                                 vjust = 2, x = x, y = y,
                                                 gp = grid::gpar(fontsize = fontsize)))
}

## applies layout matrix to grobs, then adds a white background
plot_fig <- function(grobs, layout_matrix, ...){
    p <- gridExtra::grid.arrange(grobs = grobs, layout_matrix = layout_matrix, ...)
    return(cowplot::ggdraw(p) + theme(plot.background = element_rect(fill = "white", colour = NA)))
}

## apply function to vector elements and store output as list, keeping original names
slapply <- function(v, f){
    output <- list()
    for (e in v){
        output[[e]] <- f(e)
    }
    return(output)
}

## multi map
map_one <- function(key, l_map, default=''){
    if (key %in% names(l_map)){
        return(l_map[[key]])
    }
    return(default)
}
map_multi <- function(keys, l_map, default='', auto_convert = TRUE){
    if(auto_convert & class(keys) == "factor"){keys = as.character(keys)}
    return(sapply(keys, function(k){map_one(k, l_map, default=default)}))
}

extrafont::font_import(pattern = "Arial", prompt = FALSE)
extrafont::loadfonts() ## required in chaelab-ws
theme_set(theme_void(base_family = "Arial", base_size = 8))
theme_set(theme_bw(base_family = "Arial", base_size = 8))
update_geom_defaults("text", list(size = 8/.pt))

## get all nodes descended from a given node (include_mrca tells function whether to also include mrca in output df)
## takes a data df from ggtree(tree)$data
## requires columns parent, node
all_descendants <- function(df, mrca, include_mrca = FALSE){
    if (include_mrca) {nodes <- c(mrca)} else {nodes <- c()}
    curr_parents <- c(mrca)
    while (TRUE){
        if (length(curr_parents) == 0){break}
        v.nodes <- df %>%
            dplyr::filter(parent %in% curr_parents) %>%
            dplyr::pull(node)
        if (length(v.nodes) == 0){break}
        if (all(v.nodes %in% nodes)){break}
        nodes <- c(nodes, v.nodes)
        curr_parents <- v.nodes
    }
    return(df %>% dplyr::filter(node %in% nodes))
}

## save base R plots after they've already been plotted
dev_copy <- function(fname, ...){
    for (device in list(...)){
        dev.copy(
            eval(parse(text = device)),
            paste(fname, device, sep = '.')
        )
        dev.off()
    }
}

## invert a named iterable, outputs a list of vectors
invert_named_iter <- function(input){
    l.output <- list()
    for (name in names(input)){
        value <- input[[name]]
        if (is.null(value)) { next }
        if (value %in% names(l.output)){
            l.output[[value]] <- c(l.output[[value]], name)
        } else {
            l.output[[value]] <- c(name)
        }
    }
    return(l.output)
}
