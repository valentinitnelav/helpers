make_GeodesicBuffer <- function(XY.dg, dg.step=5, dst.m, crs){
    #######################################################################################################
    ## Function to make a circle-buffer around given points (long-lat coordinates)
    ## Is different from rgeos::gBuffer() by the fact that it allows the user to create 
    ## a geodesic buffer with a width expressed in metric units.
    ## Otherwise the user needs to project and apply an Euclidian buffer with rgeos::gBuffer(),
    ## which will introduce distortions that vary greatly with latitude and the radius of the circle buffer. 
    ## See also the question addressed here:
    ## https://gis.stackexchange.com/questions/250389/euclidean-and-geodesic-buffering-in-r
    ##
    ## ARGUMENTS
    ## - XY.dg:
    ##       matrix or spatial point object (SpatialPointsDataFrame) with long-lat coordinates
    ##       (column 1 = longitude, column 2 = latitude)
    ## - dg.step:
    ##       bearing (direction) step in degrees; dictates the point density of the circle-buffer edge
    ## - dst.m:
    ##       distance (radious of circle buffer) in meters
    ## - crs:
    ##       string of class CRS-class (CRS-class {sp})
    ##       e.g. "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    ##
    ## RETURNS 
    ## SpatialPolygons as circle buffers
    #######################################################################################################
    
    # ---------------------------------------
    # Check for valid input and packages
    # ---------------------------------------
    # Check for packages
    if (!requireNamespace("geosphere", quietly = TRUE)) 
        stop("Please install package {geosphere}")
    if (!requireNamespace("data.table", quietly = TRUE)) 
        stop("Please install package {data.table}")
    if (!requireNamespace("sp", quietly = TRUE)) 
        stop("Please install package {sp}")
    
    # Check if XY.dg argument is of expected class
    is.class.ok <- inherits(XY.dg, c("SpatialPoints","SpatialPointsDataFrame","matrix"))
    if (!is.class.ok) stop("For XY.dg argument expecting 'SpatialPoints', 'SpatialPointsDataFrame' or 'matrix' class")
    # depending on class of XY.dg, get number of points
    if( inherits(XY.dg, c("SpatialPoints","SpatialPointsDataFrame")) ){
        N.points <- length(XY.dg) # if is an sp object, take length
        # Also if spatial object is projected then stop with error
        if (is.projected(XY.dg)) stop("Spatial object XY.dg is projected; expects unprojected coordinates (long-lat)")
    } else {
        # if matrix - some routine checking of coordinates values
        if (dim(XY.dg)[2] != 2) stop("Expecting XY.dg to be a 2 columns matrix (column1=longitude, column2=latitude)")
        if ( ! (range(XY.dg[,1]) %between% c(-180, 180) && range(XY.dg[,2]) %between% c(-90, 90)) )
            stop("Expects unprojected coordinates in XY.dg (longitude between -180 & 180, latitude between -90 & 90)")
        N.points <- dim(XY.dg)[1] # if is a matrix, take number of rows
    }
    
    # ===========================================================
    # A) Construct buffers as points at given distance and bearing
    # ===========================================================
    # a vector of bearings (fallows a circle)
    dg <- seq(from = 0, to = 360, by = dg.step)
    # # Construct equidistant points defining circle shapes (the "buffer points")
    # Inspired from section "Point at distance and bearing" from 
    # Robert J. Hijmans in Introduction to the ”geosphere” package at:
    # https://cran.r-project.org/web/packages/geosphere/vignettes/geosphere.pdf
    buff.XY <- geosphere::destPoint(p = XY.dg, 
                                    b = rep(dg, each = N.points), 
                                    d = dst.m)
    
    # ===========================================================
    # B) Make SpatialPolygon from the points above
    # ===========================================================
    buff.XY <- data.table(buff.XY)
    # group (split) "buffer points" by id
    # add column which indicates to which point ID from N.points each circle-buffer point belongs to
    buff.XY[, id := rep(1:N.points, times = length(dg))]
    # buff.XY[, bearing := rep(dg, each = dim(XY.dg)[1])] # for debugging reasons
    # split in a list of data tables with two columns: lon and lat
    lst <- split(buff.XY[,.(lon,lat,id)], by = "id", keep.by = FALSE)
    
    # Make SpatialPolygons out of the list of coordinates
    poly   <- lapply(lst, sp::Polygon, hole = FALSE)
    polys  <- lapply(list(poly), sp::Polygons, ID = NA)
    spolys <- sp::SpatialPolygons(Srl = polys, proj4string = CRS(crs))
    # Disaggregate (split in unique polygons)
    spolys.buff <- sp::disaggregate(spolys)
    
    return(spolys.buff)
}

## EXAMPLE
##
# require(rgeos)
# require(sp)
# require(plotKML)
# require(data.table)
# 
# # Generate a random grid-points for a (almost) global bounding box
# b.box <- as(raster::extent(120, -120, -60, 60), "SpatialPolygons")
# proj4string(b.box) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# set.seed(2017)
# pts <- sp::spsample(b.box, n=100, type="regular")
# buf1000km.geodesic <- make_GeodesicBuffer(XY.dg=pts, 
#                                           dg.step=5, 
#                                           dst.m=10^6, 
#                                           crs=CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84")))
# plot(buf1000km.geodesic)
# # save as KML 
# plotKML::kml(buf1000km.geodesic, file.name = "buf1000km.geodesic.kml")
