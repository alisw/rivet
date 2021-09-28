import yoda
def patch(path, ao):
    # remove extrapolated and average bins
    if "CLEO_2000_I526554" in path :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for i in range(0,len(ao.points())) :
            if ao.points()[i].x()==0.96 : continue
            if "d01" in path  :
                if ao.points()[i].x()==0.75 : continue
            else :
                if ao.points()[i].x()==0.72 : continue
            newAO.addPoint(ao.points()[i])
        ao=newAO
    return ao
