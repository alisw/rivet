def patch(path, ao):
    # fix bin widths
    if "ARGUS_1992_I319102" in path :
        if "d01" in path :
            step=0.1
        elif "d02" in path or "d03" in path :
            step = 1.
        elif  "d04" in path or "d05" in path or "d06" in path:
            step=0.5
        else :
            return ao
        for p in ao.points():
            p.setXErrs(step)
    return ao
