import yoda
def patch(path, ao):
    # fix hist, really a 2D histo with two equivalent x axes not a 3d one
    if ("OPAL_1993_I342766" in path) :
        if("d01" in path ) :
            points = ao.points()
            newHist=yoda.Scatter2D()
            newHist.setPath(ao.path())
            for i in range(0,len(points)) :
                x = points[i].y()
                xErrs = points[i].yErrs()
                y     = points[i].z()
                yErrs = points[i].zErrs()
                newHist.addPoint(x,y,xErrs,yErrs)    
                ao=newHist
        else :
            for p in ao.points() :
                p.setXErrs(0.2)
    return ao
