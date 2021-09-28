def patch(path, ao):
    if "OPAL_1992_I321190" in path and ("d01" in path or "d05" in path or "d06" in path or
                                        "d07" in path or "d08" in path):
        step=0.5
        if "d01" in path : step = 1.
        for p in ao.points() :
            p.setXErrs(step)
    return ao

