import math
def patch(path, ao):
    # set bin widths
    if "ARGUS_1989_I278932" in path and ( "d01" in path or "d02" in path) :
        for p in ao.points(): p.setXErrs(0.5)
    return ao
