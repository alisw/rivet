import yoda
def patch(path, ao):
    # sign issue with bin widths
    if "OPAL_1998_S3749908" in path and ("d03" in path or "d05" in path or "d07" in path or
                                         "d09" in path or "d11" in path or "d13" in path or
                                         "d15" in path):
        for p in ao.points() :
            errs=p.xErrs()
            p.setXErrs((-errs[1],-errs[0]))
    return ao
