import yoda
# removal of average bins and add bin widths, and divide by bin width
def patch(path, ao):
    if "ARGUS_1994_I354224/d01-x01-y01" in path :
        newAO = yoda.core.Scatter2D()
        newAO.setPath(ao.path())
        for i in range(0,len(ao.points())) :
            width = 0.05
            if i==19 :
                width=0.1
            elif i==20 :
                width=0.25
            # bin width
            ao.points()[i].setXErrs(width)
            # normalisation
            ao.points()[i].setY(ao.points()[i].y()/2./width)
            yerrs = ao.points()[i].yErrs()
            ao.points()[i].setYErrs((yerrs[0]/2./width,yerrs[1]/2./width))
            # remove average bin
            if i != 12 : newAO.addPoint(ao.points()[i])
        ao=newAO

    return ao
