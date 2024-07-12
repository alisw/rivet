import yoda,math
def patch(path, ao):
    if "ATLAS_2018_I1615866" in path and "d01" in path :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        newAO.setTitle(ao.title())
        for p in ao.points() :
            newPoint = yoda.Point2D(13000,p.x())
            newPoint.setXErrs(0.5)
            for (key,value) in p.errMap().items() :
                newPoint.setYErrs(value,key)
            newAO.addPoint(newPoint)
        return newAO
    else :
        return ao
