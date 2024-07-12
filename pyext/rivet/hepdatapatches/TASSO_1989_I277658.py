import yoda
def patch(path, ao):
    # bin widths (still more things need fixing)
    if "TASSO_1989_I277658" in path :
        if "d05" in path :
            for p in ao.points() :
                p.setXErrs(1)
    return ao
