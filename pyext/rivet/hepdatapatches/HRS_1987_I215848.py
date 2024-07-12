import math
def patch(path, ao):
    if "HRS_1987_I215848" in path:
        # add dummy bin width for plotting
        if "d07" in path or "d08" in path:
            for p in ao.points():
                p.setXErrs(0.5)
        # bins defined in momentum but plots in energy fraction
        # with only average given, convert momentum limits to energy fraction
        elif "d02" in path or "d03" in path or "d04" in path:
            bins=[0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9]
            ix=0
            mass=0.
            if   "d02" in path :
                mass =0.13957
            elif "d03" in path :
                mass=0.493677
            elif "d04" in path :
                mass=0.938272
            for p in ao.points():
                xmin=2.*math.sqrt(bins[ix  ]**2+mass**2)/29.
                xmax=2.*math.sqrt(bins[ix+1]**2+mass**2)/29.
                x   = 0.5*(xmin+xmax)
                wid = 0.5*(xmax-xmin)
                p.setX    ( x )
                p.setXErrs(wid)
                ix+=1
    return ao
