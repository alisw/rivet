import math
def patch(path, ao):
    # set bin widths
    if "ARGUS_1989_I268577" in path :
        if( "d01" in path) :
            for p in ao.points(): p.setXErrs(0.5)
        elif "d02" in path :
            xlim=[0,0.5,0.6,0.7,0.85,1.]
            for i in range(0,len(ao.points())) :
                ao.points()[i].setXErrs((ao.points()[i].x()-xlim[i],xlim[i+1]-ao.points()[i].x()))
    return ao
