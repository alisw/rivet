def patch(path, ao):
    # fix x values so not using bin number
    if "BESIII_2015_I1391138" in path:
        step=0.05
        if "d02" in path : step=0.1
        x = step
        for i in range(0,len(ao.points())):
            ao.points()[i].setX(x)
            ao.points()[i].setXErrs(step)
            x+=2.*step
    return ao

