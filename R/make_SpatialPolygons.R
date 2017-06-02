make_SpatialPolygons <- function(coords.lst, coords.crs)
{
    # ARGUMENTS
    # - coords.lst: 
    #       List of 2-column numeric matrices / data frames / data tables with coordinates 
    #       (column 1 = longitude, column 2 = latitude);
    #       as in sp:: Polygon, for each matrix/data frame/data table first point (row) 
    #       should equal last coordinates (row); 
    #       the status of the polygon as a hole or an island will be taken from the ring direction, 
    #       with clockwise meaning island, and counter-clockwise meaning hole.
    # - coords.crs:
    #       projection string of class CRS-class (CRS-class {sp})
    # RETURNS SpatialPolygons
    
    # Make Polygon objects out of the list of coordinates
    poly  <- lapply(coords.lst, sp::Polygon, hole = FALSE)
    # Make Polygons (note the S)
    # polys <- Map(sp::Polygons, list(poly), ID = names(poly))
    polys <- lapply(list(poly), sp::Polygons, ID = NA)
    
    # Make SpatialPolygons
    spolys <- sp::SpatialPolygons(Srl = polys, 
                                  proj4string = CRS(coords.crs))
    # Disaggregate (split in unique polygons)
    spolys <- sp::disaggregate(spolys)
    return(spolys)
    # Also check out rasterize polygons in ?rasterize for poly creation (?)
}