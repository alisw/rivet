import yoda
def patch(path, ao):
    # fix bin widths, bin ordering in hepdata gives -ve value
    if "L3_1997_I427107" in path:
        if "d06" in path :
            bins = [1.20,1.43,1.72,2.12,2.53,2.81,3.22,3.91]
        elif "d08" in path :
            bins = [1.43,1.72,2.12,2.53,2.81,3.22,3.91]
        elif "d10" in path :
            bins = [1.43,1.83,2.12,2.53,3.22,4.61]
        else :
            return ao
        for i in range(0,len(ao.points())) :
            low = ao.points()[i].x()-bins[i]
            upp = bins[i+1]-ao.points()[i].x()
            ao.points()[i].setXErrs((low,upp))
    return ao
