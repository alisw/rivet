import yoda
def patch(path, ao):
    # fix hist, really a 2D histo with two equivalent x axes not a 3d one
    if ("JADE_1990_I282847" in path and
        ("d03" in path or "d04" in path or "d05" in path)) : 
        sqrtS=35.
        if("d04" in path or "d06" in path) : sqrtS=44.
        points = ao.points()
        newHist=yoda.Scatter2D()
        newHist.setPath(ao.path())
        for i in range(0,len(points)) :
            x = 2.*points[i].x()/sqrtS
            xErrs = points[i].xErrs()
            y     = points[i].z()
            yErrs = points[i].zErrs()
            newHist.addPoint(x,y,(2.*xErrs[0]/sqrtS,2.*xErrs[1]/sqrtS),
                             (yErrs[0],yErrs[1]))    
        ao=newHist
    return ao
