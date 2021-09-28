import yoda
def patch(path, ao):
    # fix bin heights, need dividing by width
    if "OPAL_2001_I536266" in path and ("d01" in path or "d02" in path) :
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
