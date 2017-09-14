SpatialJoin_closest_within_buff_dist <- function(XY.unq, polys, dst.m, dg.step = 10){
    #######################################################################################################
    ## ! Use with care, not properly tested - gives errors for bigger than 10-20 km buffers on global scale.
    ## Function to find the closest polygon(s) to given unprojected coordinates (WGS84) within buffer distance (meters).
    ## Extracts the attribute table data corresponding to the found polygons.
    ## ARGUMENTS
    ##  XY.unq  - data.table with two columns X and Y (longitude & latitude, unique pairs)
    ##  polys   - SpatialPolygonsDataFrame object. This is used to detect the closest polygons within given buffer distance
    ##  dst.m   - buffer distance in meters
    ##  dg.step - circle point density of buffer (e.g. a vertex at each 10 dg, the default)
    ##
    ## RETURNS a data.table object with data extracted from the closest polygon.
    ## NA-s are returned for those points where no match could not be found within the given buffer.
    #######################################################################################################
    
    setnames(XY.unq, c("X","Y"))
    # Note: working only on unique pairs of coordinates can improve speed
    if( dim(unique(XY.unq, by = c("X","Y")))[1] != dim(XY.unq)[1] )
        stop("Please porvide data.table with unique pairs of coordinates; \n",
             "e.g. XY.unq <- unique(XY.dt, by = c('X','Y')); \n",
             "Very fast data.table merging can be done after using X,Y as key \n",
             "e.g. merge(x = XY.dt, y = XY.unq, by = c('X','Y'), sort = FALSE)")
    
    # create geodesic buffer arround points
    # make_GeodesicBuffer() function must be loaded
    XY.unq.buff <- make_GeodesicBuffer(XY.dg = as.matrix(XY.unq),
                                       dg.step = dg.step, # edge density of buffer (a vertex at each dg.step dg)
                                       dst.m   = dst.m,   # buffer radious
                                       crs     = proj4string(polys))
    # Intersect polys with the buffers
    # returns a list - for each polygon in polys returns the indices of the buffers intrsected
    # Order of arguments matters a lot for speed (gIntersects can be an expensive operation)
    polys.intersects.buff <- rgeos::gIntersects(spgeom1 = polys,
                                                spgeom2 = XY.unq.buff,
                                                byid = TRUE,
                                                returnDense = FALSE)
    # Note that the naming of elements in the list starts from 0 
    # (this is C or C++ way of indexing used in gdal library)
    # R indexes starting with 1, pay huge attention to such things.
    
    # Prepare the intersection results (from list to data.tble, eliminate NULL/NA cases)
    # delete original names of each element in list to avoid confusion when applying rbindlist
    names(polys.intersects.buff) <- NULL 
    # all NULL elements need to become NA (otehrwise, rbindlist does not work properly with NULL data tables)
    polys.intersects.buff[sapply(polys.intersects.buff, is.null)] <- NA
    # transform each element of the list into a data.table (setDT does not work here)
    polys.intersects.buff <- lapply(polys.intersects.buff, data.table)
    # bind tables and keep the id-s, which relates back to polys id-s (identical with point id-s)
    polys.intersects.buff <- rbindlist(polys.intersects.buff, idcol=TRUE)
    # rename columns
    setnames(polys.intersects.buff, c("poly.id", "buffORpoint.id"))
    # keep only buffer records that intersected polys
    polys.intersects.buff <- polys.intersects.buff[!is.na(buffORpoint.id)]
    
    # add an index column for polys (merging reasons)
    polys@data$poly.id <- 1:dim(polys@data)[1]
    # merge (left join) intersection results (pairs of indices) with polys data
    polys.intersects.buff <- merge(x = polys.intersects.buff,
                                   y = polys@data,
                                   by = "poly.id",
                                   all.x = TRUE)
    
    # Detect buffer id-s that intersected more than one polygons from polys
    # For such cases, it needs to detect which polygon (from those intersected) is the closest
    # and use polys data from the closest ones.
    bf.ids <- polys.intersects.buff[, .N, by = buffORpoint.id][N>=2, buffORpoint.id]
    if (length(bf.ids) != 0) {
        # select those buffers that intersected more than one polygon from polys
        buff2crop <- XY.unq.buff[bf.ids,]
        # Crop polys with buffers for cases selected above.
        # Croping a large and dense polys object reduces from computation time when computing dist2Line
        # because it reduces the number of total vertices to test for min distance.
        polys.crop <- raster::crop(x = polys, y = buff2crop) # can be very expensive operation
        # Compute shortest distance between points (center of buffers) and cropped polys;
        # also take the ID-s of closest polygons from polys (is used to subset polys)
        dist.mat <- geosphere::dist2Line(p = XY.unq[bf.ids,], line = polys.crop, distfun = distHaversine)
        # Merge buffer id-s with data from closests polys (idexed using the ID-s from dist.mat)
        closest.poly.dt <- data.table(bf.ids, polys.crop@data[dist.mat[,"ID"],])
        # move poly.id column in first position (for merging reasons)
        setcolorder(closest.poly.dt, c("poly.id", setdiff(names(closest.poly.dt), "poly.id")))
        # rename bf.ids to buffORpoint.id (for merging reasons)
        setnames(closest.poly.dt, old = "bf.ids", new = "buffORpoint.id")
        # Replace those records where a buffer intersects more than one poly in polys with
        # the closest intersected polygon data.
        # First, eliminate the respective records, then bind with the corrected ones (closests poly data prepared above)
        polys.intersects.buff <- polys.intersects.buff[!(buffORpoint.id %in% bf.ids)] # eliminate
        polys.intersects.buff <- rbindlist(list(polys.intersects.buff, closest.poly.dt)) # bind
    }
    # Prepare results
    # merge table of point XY.unq with intersection results
    XY.unq[, buffORpoint.id := 1:.N]
    XY.unq <- merge(x = XY.unq,
                    y = polys.intersects.buff,
                    by = "buffORpoint.id",
                    all.x = TRUE)
    XY.unq[, buffORpoint.id := NULL] # delete not needed column
    return(XY.unq)
}