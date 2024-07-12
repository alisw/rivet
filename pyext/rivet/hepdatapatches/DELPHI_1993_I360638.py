def patch(path, ao):
    # fix bin widths
    if "DELPHI_1993_I360638" in path and ("d02" in path or "d05" in path or "d06" in path):
        step=0.5
        for p in ao.points():
            p.setXErrs(step)
    return ao

