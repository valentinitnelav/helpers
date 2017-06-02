make_GeodesicBuffer <- function(XY.dg, dg.step, dst.m, crs)
{
    # ARGUMENTS
    # - XY.dg:
    #       matrix or spatial point object (SpatialPointsDataFrame or SpatialPoints) with long-lat coordinates
    #       (column 1 = longitude, column 2 = latitude)
    # - dg.step:
    #       bearing (direction) step in degrees; dictates the point density of the circle edge
    # - dst.m:
    #       distance (radious of circle buffer) in meters
    # - crs:
    #       projection string of class CRS-class (CRS-class {sp})
    # RETURNS SpatialPolygons circle buffer
    
    # Construct circle bearings
    dg <- seq(from = 0, to = 360, by = dg.step)
    
    # Construct circle points arround each given point
    # "Point at distance and bearing"
    # https://cran.r-project.org/web/packages/geosphere/vignettes/geosphere.pdf
    buff.XY <- geosphere::destPoint(p = XY.dg, 
                                    b = rep(dg, each = dim(XY.dg)[1]), 
                                    d = dst.m)
    buff.XY <- data.table(buff.XY)
    # add column which indicates to which point ID each circle point belongs to
    buff.XY[, id := rep(1:dim(XY.dg)[1], times = length(dg))]
    # buff.XY[, bearing := rep(dg, each = dim(XY.dg)[1])] # for debuging reasons
    setkey(buff.XY, id) # for sorting by id and improving speed (if binary search needed)
    # group (split) circle points by id
    # will result in a list of data tables with two columns: lon and lat
    lst <- split(buff.XY[,.(lon,lat,id)], by = "id", keep.by = FALSE)
    
    # From circle points to spatial polygons
    buff <- make_SpatialPolygons(coords.lst = lst,
                                 coords.crs = crs)
    return(buff)
}