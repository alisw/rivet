import math
def patch(path, ao):
    # fix bin widths, assume just transformed from the xi distribution
    if path == "/REF/DELPHI_1995_I377487/d08-x01-y01":
        bins=[0.0,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,
              3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6]
        i=len(bins)
        for p in ao.points():
            i-=1
            low = p.x()-math.exp(-bins[i])
            upp = math.exp(-bins[i-1])-p.x()
            p.setXErrs((low,upp))
    return ao

