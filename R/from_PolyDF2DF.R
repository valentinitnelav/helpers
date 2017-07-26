from_PolyDF2DF <- function(my.polyDF)
{
    # Or left join results of broom::tidy() with my.polyDF@data  by "id" column 
    # example with ggplot2::fortify() at https://github.com/Robinlovelace/Creating-maps-in-R
    
    # ----------------------------------------
    # Function to convert from SpatialPolygonsDataFrame to data frame while keeping attributes
    # my.polyDF = SpatialPolygonsDataFrame object
    # returns a data frame object
    # reference: https://gis.stackexchange.com/questions/169599/how-to-extract-all-the-polygon-coordinates-from-a-spatialpolygonsdataframe
    # ----------------------------------------
    
    # transform my.polyDF to SpatialLinesDataFrame object
    lines <- as(my.polyDF, "SpatialLinesDataFrame")
    # transform the lines to data frame - preserves data associated with the lines 
    my.df <- as.data.frame(as(lines, "SpatialPointsDataFrame"))
    # Transform my.polyDF directly to data frame, but this does not preserves data associated with each poly
    df.tidy <- broom::tidy(my.polyDF)
    # Merge the two data frames to have both the coordinates and the attributes associated with them
    results <- cbind(df.tidy, my.df[,!names(my.df) %in% c("coords.x1",  "coords.x2")])
    return(results)
}