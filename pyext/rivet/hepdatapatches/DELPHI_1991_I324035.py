def patch(path, ao):
    # fix bin widths
    if "DELPHI_1991_I324035" in path:
        step=0.5
        if("d05" in path) : step=1.
        for p in ao.points():
            p.setXErrs(step)
    return ao

