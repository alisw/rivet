def patch(path, ao):
    # fix bin widths
    if "ARGUS_1993_S2653028" in path and ("d07" in path or "d08" in path or "d09" in path or
                                          "d10" in path or "d11" in path):
        for p in ao.points():
            p.setXErrs(0.5)
    return ao
