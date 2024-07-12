import yoda
def patch(path, ao):
    # fix bin widths
    if path == "/REF/L3_1995_I381046/d01-x01-y01" :
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
