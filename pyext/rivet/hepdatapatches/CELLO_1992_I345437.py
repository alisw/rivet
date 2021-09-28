
def patch(path, ao):
    if "CELLO_1992_I345437" in path :
        for p in ao.points():
            if "d01" in path  :
                if(p.x()<1.5) :
                    p.setXErrs(0.0125)
                else :
                    p.setXErrs(0.05)
            else :
                p.setXErrs(0.025)
    return ao

