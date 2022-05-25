def patch(path, ao):
    # fix bin widths
    if "ARGUS_1993_S2789213" in path and ("d01" in path or "d02" in path or "d03" in path or
                                          "d16" in path or "d17" in path):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
