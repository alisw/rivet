import yoda
def patch(path, ao):
    # bin widths
    if path == "/REF/OPAL_1997_I440103/d04-x01-y01":
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
