def patch(path, ao):
    # fix bin widths
    if "/REF/HRS_1986_I18688/d01"  in path:
        for p in ao.points() :
            if p.x()<0.4 : p.setXErrs(0.05)
            else :         p.setXErrs(0.15)
    return ao
