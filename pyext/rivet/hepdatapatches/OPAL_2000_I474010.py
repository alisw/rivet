import yoda
def patch(path, ao):
    # bin widths
    if "OPAL_2000_I474010" in path and ( "d01" in path or "d02" in path or "d03" in path) :
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
