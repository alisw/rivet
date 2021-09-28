import math
def patch(path, ao):
    if "AMY_1990_I295160" in path and ("d01" in path or "d02" in path):
        for p in ao.points():
            p.setXErrs(1)
    return ao
