def patch(path, ao):
    # fix bin widths
    if "DELPHI_1992_I334948" in path:
        step=1.
        for p in ao.points():
            p.setXErrs(step)
    return ao

