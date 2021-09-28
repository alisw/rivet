import yoda
def patch(path, ao):
    # sign issue with bin widths
    if path == "/REF/OPAL_1998_I474012/d01-x01-y01":
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
