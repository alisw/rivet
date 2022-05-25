import yoda
def patch(path, ao):
    # bin widths
    if "OPAL_1998_S3780481" in path and "d09" in path :
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
