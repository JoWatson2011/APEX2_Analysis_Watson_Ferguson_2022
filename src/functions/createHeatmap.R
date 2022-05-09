hclust_heatmap <- function(df_long, k = NULL, h = NULL) {
  set.seed(1)
  tree <- hclust(dist(df_long))
  
  coupe <- cutree(tree, k = k, h = h)
  
  coupe.or <- coupe[tree$order]
  coupe.out <- rep(NA, length(coupe))
  j <- 1 #
  k <- coupe.or[1]
  for (i in 1:length(coupe)) {
    if (coupe.or[i] == k) {
      next
    } else {
      coupe.out[which(coupe == k)] <- j
      j <- j + 1
      k <- coupe.or[i]
    }
  }
  coupe.out[is.na(coupe.out)] <- j
  names(coupe.out) <- names(coupe)
  coupe.out
}

create_heatmap <- function(df_long, den_cl) {

  set.seed(1)
  plot_heatmap <- function() gplots::heatmap.2(x = df_long,
                                               Rowv = as.dendrogram(hclust(dist(df_long))),
                                               RowSideColors=paletteer_d("ggsci::default_igv", max(den_cl))[den_cl],
                                               col=colorRampPalette(c('dodgerblue4','lightskyblue2','white','pink','darkred')),
                                               main="Clusters numbered low>high from bottom>top",
                                               trace="none",
                                               density.info="none",
                                               labRow = "",
                                               cexCol= 0.75,
                                               srtCol = 45,
                                               key.title = "Z-score",
                                               offsetCol = 0)
}
