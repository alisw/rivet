def patch(path, ao):
    if "ALEPH_1991_S2435284" in path:
        step = 1.
        if "d01" not in path : step=0.5
        for p in ao.points():
            p.setXErrs(step)
    return ao
