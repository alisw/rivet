def patch(path, ao):
    # fix bin widths
    if "CRYSTAL_BALL_1991_I297905" in path and ("d01" in path or "d02" in path):
        for p in ao.points():
            if "d01" in path :
                p.setX(10.05)
                p.setXErrs(0.65)
            else :
                p.setXErrs(0.5)
    return ao
