import math
def patch(path, ao):
    if "CELLO_1982_I12010" in path:
        step = math.pi/100.
        x = step
        if "d03" in path or "d04" in path: x+=2.*step
        for p in ao.points():
            p.setX(x)
            p.setXErrs(step)
            x+=2.*step
    return ao
