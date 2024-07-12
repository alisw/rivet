# bin widths
def patch(path, ao):
    if "L3_1990_I298078" in path and "d01" in path:
        for i in range(0,len(ao.points())) :
            if(i<6) :
                ao.points()[i].setXErrs(5e-3)
            elif(i==6) :
                ao.points()[i].setXErrs((5e-3,1e-2))
            else :
                ao.points()[i].setXErrs(0.01)
    return ao
