import yoda
def patch(path, ao):
    # fix bin heights, need dividing by width
    if path == "/REF/OPAL_2003_I599181/d01-x01-y01" :
        for p in ao.points() :
            width = p.xErrs()[0]+p.xErrs()[1]
            yErrs = p.yErrs()
            for val in yErrs:
                val /= width
            p.setY(p.y()/width)
            p.setYErrs(yErrs)
    return ao
