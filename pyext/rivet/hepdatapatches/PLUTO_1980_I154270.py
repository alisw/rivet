import yoda,math
# removal of average bins and bin number -> cms energy
def patch(path, ao):
    if path == "/REF/PLUTO_1980_I154270/d01-x01-y01" :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for i in range(0,len(ao.points())) :
            if ao.points()[i].x()!=30.75 :
                newAO.addPoint(ao.points()[i])
        ao=newAO

    return ao
