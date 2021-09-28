import yoda
def patch(path, ao):
    # sign issue with bin widths (still need to fix ordering issue)
    if path == "/REF/TASSO_1989_I267755/d05-x01-y01":
        for p in ao.points() :
            errs=p.xErrs()
            if(errs[0]<0.) :
                p.setXErrs((-errs[1],-errs[0]))
    return ao
