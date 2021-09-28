import yoda,math
# removal of average bins and bin number -> cms energy
def patch(path, ao):
    if path == "/REF/PLUTO_1979_I142517/d01-x01-y01" :
        bins = [22.,27.6,30.0,31.6]
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for i in range(0,len(ao.points())) :
            if i != 4 :
                ao.points()[i].setX(bins[i])
                ao.points()[i].setXErrs(0.)
                newAO.addPoint(ao.points()[i])
        ao=newAO

    return ao
