import yoda
def patch(path, ao):
    if "MARKII_1985_I207785" in path:
        # set bin width
        if "d01" in path:
            for p in ao.points():
                p.setXErrs(0.5)
        # really only normal hist with mapped axis defined
        elif "d02" in path or "d04" in path :
            newAO = yoda.core.Scatter2D()
            newAO.setPath(ao.path())
            points = ao.points()
            for i in range(0,len(points)) :
                x = 2.*points[i].x()/29.
                xErrs = points[i].xErrs()
                y     = points[i].z()
                yErrs = points[i].zErrs()
                newAO.addPoint(x,y,(2.*xErrs.minus/29.,2.*xErrs.plus/29.),yErrs)
            ao=newAO
    return ao
