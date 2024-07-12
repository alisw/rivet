import yoda
def patch(path, ao):
    if "MD1_1994_I362599" in path:
        if "d02" in path or "d04" in path :
            for p in ao.points() : p.setXErrs(0.5)
        elif "d03" in path :
            newHists=[]
            for i in range(1,3) :
                newHists.append(yoda.Scatter2D())
                newHists[-1].setPath(ao.path().replace("y01","y0%s"%i))
            points = ao.points()
            for i,p in enumerate(points):
                newHists[i].addPoint(p)
            return newHists
    return ao
