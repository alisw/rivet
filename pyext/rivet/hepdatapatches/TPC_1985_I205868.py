import yoda
def patch(path, ao):
    # fix type of hist 2d not 3d, 2 equivalent x axes
    if "TPC_1985_I205868" in path:
        points = ao.points()
        newHist=yoda.Scatter2D()
        newHist.setPath(ao.path())
        for i in range(0,len(points)) :
            x = 2.*points[i].x()/29.
            xErrs = points[i].xErrs()
            xErrs = (xErrs[0]*2./29.,xErrs[1]*2./29.)
            y     = points[i].z()
            yErrs = points[i].zErrs()
            newHist.addPoint(x,y,xErrs,yErrs)    
        ao=newHist
    return ao
