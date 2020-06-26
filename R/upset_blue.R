library(ComplexHeatmap)
library(SuperExactTest)
library(ggpubr)

#--------------------------------------------------------------------------------

SET_input <- function(m){
  out <- lapply(colnames(m),function(x){
    in_pop <- m[,x] == 1
    rownames(m)[in_pop]
  })
  names(out) <- colnames(m)
  out}


blue_fun <- function(x){
  ref_val <- seq(0,ceiling(max(x)))
  col <- RColorBrewer::brewer.pal(9,"Blues")
  pal <- colorRampPalette(col)(length(ref_val))
  circlize::colorRamp2(ref_val,pal)}



size_order <- function(m){  # m is a combination matrix for UpSet
  sorter <- data.frame(
    idx    = 1:length(comb_name(m)),
    code   = comb_name(m),
    size   = comb_size(m),
    degree =  comb_degree(m),
    level  =  lapply(strsplit(comb_name(m),""),
                     function(x) {
                       v <- letters[1:4][as.numeric(x) == 1]
                       paste0(v, collapse = "")}) %>% unlist()
  )
  with(sorter,idx[order(-size,degree,level)])
}


numbered_bar <- function(m, # m is a combination matrix for UpSet
                         comb_order = size_order(m),
                         fill = "black",
                         height = unit(12, "cm"),
                         n_offset =1.5,
                         col_fun = NULL
                         ){

  if(is.numeric(fill)){
    fill <- col_fun(fill)
  }

  HeatmapAnnotation(
  "Intersection Size" = anno_barplot(
    comb_size(m),
    border = FALSE,
    gp = gpar(fill = fill),
    height = height,
    axis_param = list(side = "left",
                      gp=gpar(fontsize=12)),
  ),

  annotation_name_side = "left",
  annotation_name_rot = 90)
}


add_numbers_lgd <- function(m,
  n_offset = NULL,
  negLogP = NULL,
  comb_order = size_order(m)) {
  max_y <- max(comb_size(m))
  max_x <-  dim(m)[2]
  n_offset <- max_y/25
  decorate_annotation("Intersection Size", {
  # value on x-axis is always 1:ncol(mat)
  # while values on y-axis is the value after column reordering
  grid.text(label = comb_size(m[comb_order]),
            x= 1:max_x,
            y = comb_size(m[comb_order]) + n_offset,
            just = "top",
            default.units = "native",
            draw= TRUE)
  lgd <-  Legend(col_fun   = blue_fun(negLogP),
                   title     = expression( -"Log"[10] ( italic(P) )),
                   direction = "horizontal")
  draw(lgd,
       x = unit(max_x,"native"),
       y = unit(max_y,"native"),
       just = c("right", "top"))
  })
}

set_size_text <- function(m, fontsize = 12){
  HeatmapAnnotation(
    size =  anno_text( set_size(m), gp = gpar(fontsize = fontsize)),
    which = "row")
}


plot_blue <- function(m, pvalues, ...){
  zeros <- grep(0,pvalues)
  pvalues[zeros] <- 1e-320

  comb <- make_comb_mat(m, mode = "intersect")
  comb <- comb[comb_degree(comb) >1]

  comb_m <- comb[size_order(comb)]
  negLogP <- -log10(pvalues[comb_name(comb_m)])

  UpSet(
    comb_m,
    comb_order = size_order(comb_m),
    set_order = pops,
    pt_size = unit(20, "points"),
    lwd = 5,
    height = unit(0.8, "cm")*nrow(comb_m),
    #width  = unit(0.8, "cm")*ncol(comb_m),
    top_annotation = numbered_bar(
      comb_m,
      height = unit(12, "cm"),
      n_offset = 1.5,
      fill = negLogP,
      col_fun = blue_fun(negLogP)),
    right_annotation = set_size_text(comb_m)
  ) %>% draw()

  add_numbers_lgd(comb_m, negLogP = negLogP)
}

