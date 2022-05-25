import yoda
def patch(path, ao):
    # bin widths
    if "OPAL_1994_S2927284" in path and "d04" in path:
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
