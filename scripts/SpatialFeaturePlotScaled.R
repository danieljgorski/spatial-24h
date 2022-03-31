# This function adapts the SpatialFeaturePlot() function, but corrects so that,
# features with different scales on two images are plotted with the same scale.
# this is specially adapted for spatial-transcriptomics-24h-IR using "Sham",
# and "IR" images.

SpatialFeaturePlotScaled <- function(object,
                                     group,
                                     group.1,
                                     group.2,
                                     feature_of_interest,
                                     from.meta.data = TRUE,
                                     group.1.title,
                                     group.2.title,
                                     legend.title){
  if (from.meta.data == TRUE) {
    
    #subset based on input group
    Idents(object) <- group
    object_sub1 <- subset(object, idents=group.1)
    object_sub2 <- subset(object, idents=group.2)
    
    # find max range for feature of interest in each subgroup
    sub1_max <-max(object_sub1[[feature_of_interest]])
    sub2_max <-max(object_sub2[[feature_of_interest]])
    
    # set colors
    replicants <- c("#5E4FA2",
                    "#388FBA", 
                    "#A7DBA4",
                    "#E0F298",
                    "#F5FBB0",
                    "#FDD27F",
                    "#F78850",
                    "#D8434D",
                    "#AB1044")
    
    # plot
    a <- SpatialFeaturePlot(object, 
                            features = feature_of_interest,
                            images = group.1,
                            pt.size.factor = 1.8) + 
      ggtitle(group.1.title) +
      theme(legend.position = "right", 
            title = element_text(size = 16)) & 
      scale_fill_gradientn(colours=replicants,
                           name=legend.title,
                           limits = c(0, max(c(sub1_max, sub2_max))))
    b <- SpatialFeaturePlot(object,
                            features = feature_of_interest,
                            images = group.2,
                            pt.size.factor = 1.8,) + 
      ggtitle(group.2.title) +
      theme(legend.position = "right", 
            title = element_text(size = 16)) & 
      scale_fill_gradientn(colours=replicants,
                           name=legend.title,
                           limits = c(0, max(c(sub1_max, sub2_max))))
    print(a + b + plot_layout(guides = "collect"))
    
  } else {
    
    #subset based on input group
    Idents(object) <- group
    object_sub1 <- subset(object, idents=group.1)
    object_sub2 <- subset(object, idents=group.2)
    
    # find max range for feature of interest in each subgroup
    sub1_max <-max(GetAssayData(object = object_sub1,
                                slot = "data")[feature_of_interest,])
    sub2_max <-max(GetAssayData(object = object_sub2,
                                slot = "data")[feature_of_interest,])
    
    # set colors
    replicants <- c("#5E4FA2",
                    "#388FBA", 
                    "#A7DBA4",
                    "#E0F298",
                    "#F5FBB0",
                    "#FDD27F",
                    "#F78850",
                    "#D8434D",
                    "#AB1044")
    
    # plot
    a <- SpatialFeaturePlot(object,
                            features = feature_of_interest,
                            images = group.1,
                            pt.size.factor = 1.8,
                            alpha = c(0.01,1)) + 
      ggtitle(group.1.title) +
      theme(legend.position = "right", 
            legend.title = element_text(face = "italic", size = 16),
            title = element_text(size = 16)) & 
      scale_fill_gradientn(colours = replicants,
                           name = feature_of_interest,
                           limits = c(0, max(c(sub1_max, sub2_max))))
    b <- SpatialFeaturePlot(object,
                            features = feature_of_interest,
                            images = group.2,
                            pt.size.factor = 1.8,
                            alpha = c(0.01,1)) + 
      ggtitle(group.2.title) +
      theme(legend.position = "right", 
            legend.title = element_text(face = "italic", size = 16),
            title = element_text(size = 16)) & 
      scale_fill_gradientn(colours = replicants,
                           name=feature_of_interest,
                           limits = c(0, max(c(sub1_max, sub2_max))))
    print(a + b + plot_layout(guides = "collect"))
  }
}

