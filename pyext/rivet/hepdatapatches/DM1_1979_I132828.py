def patch(path, ao):
    # fix bin values, from bin number to cms energy
    x = [982.5,1016.,1056.,1098.]
    if path == "/REF/DM1_1979_I132828/d01-x01-y01":
        for i in range(0,len(ao.points())):
            ao.points()[i].setX(x[i])
            if i==0 :
                ao.points()[i].setXErrs(19.5)
            elif i==1 :
                ao.points()[i].setXErrs(8.)
            else :
                ao.points()[i].setXErrs(0.)
    return ao
