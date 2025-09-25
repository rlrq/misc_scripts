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
    dir.create(dir, showWarnings = FALSE)
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
