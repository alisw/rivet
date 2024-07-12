def patch(path, ao):
    if "VENUS_1995_I392360" in path:
        step = 0.05
        if "d01" in path : step=0.025
        for p in ao.points():
            p.setXErrs(step)
    return ao
