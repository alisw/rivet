import yoda
def patch(path, ao):
    # no bin widths, need to add them
    if "DELPHI_1990_I297698" in path :
        for i in range(0,len(ao.points())) :
            if i<=4  :
                ao.points()[i].setXErrs(0.005)
            elif i==5 :
                ao.points()[i].setXErrs((0.005,0.01))
            elif i<=13 :
                ao.points()[i].setXErrs(0.01)
            else :
                ao.points()[i].setXErrs(0.02)
    return ao
