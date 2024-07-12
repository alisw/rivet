import yoda,math
# removal of average bins and bin number -> cms energy
def patch(path, ao):
    if path == "/REF/PLUTO_1981_I165122/d05-x01-y01" :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for i in range(0,len(ao.points())) :
            if ao.points()[i].x()!=0.53 :
                newAO.addPoint(ao.points()[i])
        ao=newAO

    return ao
