def patch(path, ao):
    # change x-axis to cms energy not bin number
    x = [1900,1920,1925,1940,1950,1960,1975,1980,2000]
    # fix bin widths
    if "CMD3_2016_I1385598" in path and "d01" in path :
        for i in range(0,len(ao.points())) :
            if(i<len(x)) : ao.points()[i].setX(x[i])
    return ao
