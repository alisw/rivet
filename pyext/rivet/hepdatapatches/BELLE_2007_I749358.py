def patch(path, ao):
    if "BELLE_2007_I749358" in path:
        step = 0.025
        if "d01" in path : step=0.0025
        for p in ao.points():
            p.setXErrs(step)
    return ao
