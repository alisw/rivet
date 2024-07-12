import yoda,math
# removal of average bins and replace bin no with x values
def patch(path, ao):
    if "BABAR_2016_I1391152" in path :
        mB = 5.27966
        mK = .89167
        bins=[(1.,6.),(0.1,2),(2,4.3),(4.3,8.12),(10.11,12.89),(14.21,(mB-mK)**2)]
        output=[yoda.core.Scatter2D(),yoda.core.Scatter2D()]
        output[0].setPath(ao.path())
        output[1].setPath(ao.path().replace("x01","x02"))
        for i, p in enumerate(ao.points()) :
            x  = 0.5*(bins[i][1]+bins[i][0])
            dx = 0.5*(bins[i][1]-bins[i][0])
            p.setX(x)
            p.setXErrs(dx)
            if i!=0 :
                output[0].addPoint(p)
            else :
                output[1].addPoint(p)
        return output
    return ao
