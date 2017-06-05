extract_CharBetween <- function(my.txt, left, right)
{
    # Function to extract characters between given characters
    # __ Argumets:
    # my.txt = Character variable / line of text from where to extract
    # left   = characters - left edge to extract
    # right  = characters - right edge to extract; 
    # note: use right = "" to flag that there is no right part
    # __ Returns value:
    # returns character
    
    # if the riht part is "" (nothing)
    if (right=="")
    {
        # split the text in two, take only the last (2nd) part
        split_1_2nd_part <- strsplit(x = my.txt, split = left)[[1]][2]
        return(split_1_2nd_part)
    } else
    {
        # split the text in two, take only the last (2nd) part
        split_1_2nd_part <- strsplit(x = my.txt, split = left)[[1]][2]
        # split the 2nd part from above in two and take only the first part
        split_2_1st_part <- strsplit(x = split_1_2nd_part, split = right)[[1]][1]
        return(split_2_1st_part)
    }
}