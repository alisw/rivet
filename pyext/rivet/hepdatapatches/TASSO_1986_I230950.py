# fix bin widths (one in hepdata are for p not x)
def patch(path, ao):
    if path == "/REF/TASSO_1986_I230950/d02-x01-y01" :
        bins = [0.7,1.0,1.5,2.0,3.0,4.0,6.0,10.0,17.0]
        i=0
        for p in ao.points() :
            low = p.x()-2.*bins[i]/34.4
            upp = 2.*bins[i+1]/34.4-p.x()
            p.setXErrs((low,upp))
            i+=1
    return ao
