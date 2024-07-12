def patch(path, ao):
    # fix bin widths
    if "DELPHI_1991_I301657" in path and "d02" in path:
        for p in ao.points():
            p.setXErrs(1.)
    return ao

