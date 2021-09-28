import yoda
def patch(path, ao):
    # fix bin heights, need dividing by width
    if "OPAL_2000_I513476" in path and "d13" in path :
        for p in ao.points() :
            p.setXErrs(1.)
    return ao
