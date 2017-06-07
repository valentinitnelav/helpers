sample_SpatialPointsWithCircles <- function(sp.polygons, 
                                            sp.points, 
                                            N.random.circles, 
                                            radius,
                                            containing,
                                            my.seed = 2017L)
{
    # ----------------------------------------
    # Function to sample given points with circles at given random locations within given polygon
    # ___ Arguments
    # sp.polygons       - Spatial polygons within sampling takes place. Must have the same CRS as sp.points.
    #                       The CRS must be a projected one (no long-lat, but in meters)
    # sp.points         - Spatial points that are sampled. Must have the same CRS as sp.polygons.
    # N.random.circles  - Number of circles used to sample within the area of sp.polygons.
    # radius            - Radius of a sampling circle in meters. Same details as for width argument in rgeos::gBuffer()
    # containing        - Exact number (integer) of sampled points to be contained by a sampling circle.
    # my.seed           - Random seed (integer).

    # ___ Returns
    # A list with a data frame and a spatial object
    # circle.sample  - data frame with points sampled and associated data
    # circles.OK     - spatial object with the circles used for sampling
    # ----------------------------------------
    
    # generate random points and make circles
    set.seed(seed = my.seed)
    pts.random <- sp::spsample(sp.polygons, 
                               n    = N.random.circles, 
                               type = "random", 
                               iter = 10)
    # make circles
    circles <- rgeos::gBuffer(spgeom = pts.random, 
                              byid   = TRUE, 
                              width  = radius)
    
    # keep only circles inside land
    # disolve polygons to a single polygon
    sp.polygons.agg <- raster::aggregate(x = sp.polygons)
    IsCircleInLand <- rgeos::gContains(spgeom1 = sp.polygons.agg,
                                       spgeom2 = circles,
                                       byid    = TRUE)
    
    # subset circles
    circles.InLand <- circles[IsCircleInLand,]
    
    # intersetc points with in-land circles
    IsPointInCircle <- rgeos::gContains(spgeom1 = circles.InLand, 
                                        spgeom2 = sp.points, 
                                        byid    = TRUE)
    
    # keep only circles that respect condition (e.g. only containing X points)
    IsCircleOK <- colSums(IsPointInCircle) == containing
    # sum(IsCircleOK)
    circles.OK <- circles.InLand[IsCircleOK,]
    
    # get polygons intersected by selected circles
    circ.over.polys <- sp::over(x = circles.OK, sp.polygons[,c(1,2)])
    circ.over.polys$circle.id <- 1:dim(circ.over.polys)[1]
    
    # get point IDs from within selected circles
    IsPointInCircleOK <- IsPointInCircle[, IsCircleOK, drop=FALSE]
    
    pts.sampled <- which(IsPointInCircleOK, arr.ind = T)
    pts.sampled <- cbind(pts.sampled, as.integer(rownames(pts.sampled)))
    colnames(pts.sampled)[2:3] <- c("circle.id", "point.ID")
    
    # prepare final results
    circle.sample <- merge(x  = pts.sampled, 
                           y  = circ.over.polys, 
                           by = "circle.id")
    circle.sample[,"row"] <- NULL
    colnames(circle.sample)[3:4] <- c("geoentity.ID", "geoentity.name")
    
    results <- list(circle.sample, circles.OK)
    names(results) <- c("circle.sample", "circles.OK")
    return(results)
}