import yoda,math
def patch(path, ao):
    # weird data in hepdata review of identified particles but not in the paper entry!
    if "TASSO_1984_I194774" in path and "d01" in path:
        aopath=ao.path()
        ao=yoda.Scatter2D()
        ao.setPath(aopath)
        # fix bin widths
        x   = [0.375,0.525,0.675,0.825]
        y   = [10.3,3.7,7.2,2.1]
        err = [5.3,3.4,2.5,1.5]
        # PDG 2020
        br  = 0.045
        br_err = 0.004
        step=0.075
        for i in range(0,len(x)) :
            ao.addPoint(x[i],y[i]/br,step,y[i]/br*math.sqrt((err[i]/y[i])**2+(br_err/br)**2))
    return ao
