
import yoda
def patch(path, ao):
    # last bin not properly normalised (strip it off)
    if "MARKII_1988_I246184" in path and ("d11" in path or "d12" in path or "d14" in path or
                                          "d29" in path or "d30" in path or "d32" in path or
                                          "d47" in path or "d48" in path or "d50" in path ) :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        points = ao.points()
        for i in range(0,len(points)-1) :
            newAO.addPoint(points[i])
        ao=newAO
    return ao
    
