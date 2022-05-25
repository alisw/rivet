import math
def patch(path, ao):
    # set bin widths
    if "ARGUS_1989_I262551" in path and ("d03" in path or "d04" in path or
                                         "d05" in path or "d06" in path):
        for p in ao.points() :
            p.setXErrs(0.5)
    return ao
